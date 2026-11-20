import os
import random
import shutil
import subprocess
import sys
from itertools import product
from pathlib import Path
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import plotly.graph_objs as go

def blastn(db_path, work_folder, blastn_exe='blastn', n_cores=5, n_seqs=1, subset_size=100):
    db_folder = db_path / 'db'
    db_taxonomy = db_path / 'db_taxonomy.parquet.snappy'
    db_fasta = work_folder / 'tmp_db.fasta'
    query_fasta = work_folder / 'tmp_query.fasta'
    subset_folder = work_folder / 'tmp_subsets'
    db_blastn_csv = work_folder / 'tmp_db_blastn.csv'

    # Skip if BLAST was already run completely
    if db_blastn_csv.exists():
        print(f'{db_blastn_csv.name} already exists - skipping BLAST')
        return db_blastn_csv

    # 1) Create the query fasta (skipped if it already exists from a previous start)
    if not query_fasta.exists():
        # Create fasta from database
        subprocess.call(['blastdbcmd', '-db', str(db_folder), '-entry', 'all', '-outfmt', '%f', '-out', str(db_fasta)])
        if not db_fasta.exists():
            print('Could not export the database to fasta')
            return

        # Limit to n_seqs random sequences per species (only proper species names, no "sp.")
        taxonomy_df = pd.read_parquet(db_taxonomy).fillna('')
        taxonomy_df = taxonomy_df[taxonomy_df['species'] != '']
        taxonomy_df = taxonomy_df[~taxonomy_df['species'].str.contains(r'\bsp\.', regex=True)]
        selected = set(taxonomy_df.sample(frac=1, random_state=1).groupby('species').head(n_seqs)['Accession'])

        records = [record for record in SeqIO.parse(db_fasta, 'fasta') if record.id in selected]
        db_fasta.unlink()
        if len(records) == 0:
            print('No query sequences found - do the fasta headers match the "Accession" column?')
            return
        SeqIO.write(records, query_fasta, 'fasta')
        print(f'{len(records):,} query sequences ({n_seqs} per species)')

    # 2) Split the query fasta into subsets (skipped if they already exist)
    if not subset_folder.exists():
        os.makedirs(subset_folder)
        records = list(SeqIO.parse(query_fasta, 'fasta'))
        for i in range(0, len(records), subset_size):
            subset_fasta = subset_folder / f'subset_{i // subset_size + 1}.fasta'
            SeqIO.write(records[i:i + subset_size], subset_fasta, 'fasta')

    # 3) BLAST each subset against the full database; finished subsets are skipped
    subset_fastas = sorted(subset_folder.glob('subset_*.fasta'), key=lambda f: int(f.stem.split('_')[1]))
    for subset_fasta in tqdm(subset_fastas, desc='BLAST subsets'):
        subset_csv = subset_fasta.with_suffix('.csv')
        if subset_csv.exists():
            continue

        # Write to a temporary file first, so an interrupted subset is never mistaken for a finished one
        subset_tmp = subset_fasta.with_suffix('.running')
        result = subprocess.call([blastn_exe, '-task', 'megablast', '-db', str(db_folder), '-query', str(subset_fasta),
                                  '-num_threads', str(n_cores), '-max_target_seqs', '250',
                                  '-outfmt', '6 delim=;; qseqid sseqid pident evalue qcovs', '-out', str(subset_tmp)])
        if result != 0:
            print(f'BLAST failed for {subset_fasta.name} - start again to continue from here')
            if subset_tmp.exists():
                subset_tmp.unlink()
            return
        subset_tmp.rename(subset_csv)

    # 4) Combine all subset results into one file and remove the subsets
    with open(db_blastn_csv, 'w') as f_out:
        for subset_fasta in subset_fastas:
            with open(subset_fasta.with_suffix('.csv')) as f_in:
                f_out.write(f_in.read())
    shutil.rmtree(subset_folder)

    return db_blastn_csv

def optimal_threshold(positives, negatives, fp_weight=1, step=0.1):
    """
    Threshold that best separates two similarity distributions of one taxonomic level (barcoding gap):
    - positives: similarity to the closest hit of the SAME taxon (should be assigned to this taxon)
    - negatives: similarity to the closest hit of ANOTHER taxon (would be a wrong assignment)
    False negatives (FN) = positives below the threshold, false positives (FP) = negatives at or above it.
    The threshold with the lowest weighted error (fp_weight * FP + FN) is chosen. If several thresholds are
    equally good (e.g. a gap between both distributions), the middle one is used.
    """
    candidates = np.round(np.arange(50, 100 + step, step), 1)
    curve = []
    for t in candidates:
        fn = (positives < t).mean() * 100
        fp = (negatives >= t).mean() * 100
        curve.append([t, fn, fp, fp_weight * fp + fn])
    curve = pd.DataFrame(curve, columns=['Threshold', 'False negatives (%)', 'False positives (%)', 'Weighted error'])
    best = curve[curve['Weighted error'] == curve['Weighted error'].min()]
    optimum = best.iloc[len(best) // 2]
    return optimum['Threshold'], optimum['False negatives (%)'], optimum['False positives (%)'], curve

def summarise_levels(res_df, test_levels, fp_weight=1, fallback=True):
    # Summary and optimal threshold per level
    summary = []
    curves = []
    for level in test_levels:
        # Other taxon (negatives): only queries with relatives and below 100% similarity
        values = res_df[level].dropna()
        has_relatives = res_df.loc[values.index, f'{level} (has relatives)'] == True
        identical = values == 100
        negatives = values[has_relatives & ~identical]
        if len(negatives) == 0 and fallback == True:
            negatives = values

        # Same taxon (positives)
        positives = res_df[f'{level} (same taxon)'].dropna()

        if len(positives) > 0 and len(negatives) > 0:
            threshold, fn, fp, curve = optimal_threshold(positives.values, negatives.values, fp_weight)
            curve.insert(0, 'Level', level)
            curves.append(curve)
        else:
            threshold, fn, fp = np.nan, np.nan, np.nan

        summary.append([level, len(values),
                        round(identical.mean() * 100, 1),         # % identical to another taxon
                        round((~has_relatives).mean() * 100, 1),  # % isolated lineages
                        len(positives), len(negatives),
                        threshold, round(fn, 1), round(fp, 1),
                        positives.median(), negatives.median(), negatives.quantile(0.95)])
    summary = pd.DataFrame(summary, columns=['Level', 'n queries', '% identical (100%)', '% isolated lineages',
                                             'n same taxon', 'n other taxon',
                                             'Optimal threshold', 'False negatives (%)', 'False positives (%)',
                                             'Median same taxon', 'Median other taxon', '95% other taxon'])
    return summary, curves

def leave_one_out(blast_csv, taxonomy_parquet, output_xlsx, min_coverage=80, min_relatives=2, taxon_level='class', min_group_size=10, fp_weight=1):
    tax_levels = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    test_levels = ['species', 'genus', 'family', 'order', 'class']
    parent_level = {'species': 'genus', 'genus': 'family', 'family': 'order', 'order': 'class', 'class': 'phylum'}
    child_level = {'species': None, 'genus': 'species', 'family': 'genus', 'order': 'family', 'class': 'order'}

    # BLAST hits (with query coverage if available)
    df = pd.read_csv(blast_csv, engine='python', sep=';;', header=None)
    if len(df.columns) == 5:
        df.columns = ['Query', 'Accession', 'Similarity', 'E-value', 'Coverage']
        # Remove short local alignments (they can reach 100% similarity with distant taxa)
        df = df[df['Coverage'] >= min_coverage]
    else:
        df.columns = ['Query', 'Accession', 'Similarity', 'E-value']
        print('Note: no query coverage in the BLAST results - rerun BLAST to use the coverage filter.')

    # Reference taxonomy of each hit
    taxonomy_df = pd.read_parquet(taxonomy_parquet).fillna('')
    df = df.merge(taxonomy_df[['Accession'] + tax_levels], on='Accession', how='left')

    # Lookups: accession -> species, species -> full taxonomy
    accession_species = dict(zip(taxonomy_df['Accession'], taxonomy_df['species']))
    species_taxonomy = taxonomy_df.drop_duplicates('species').set_index('species', drop=False)[tax_levels]

    # Number of members of each parent taxon in the database (e.g. species per genus)
    n_relatives = {}
    for level in test_levels:
        parent = parent_level[level]
        sub = taxonomy_df[(taxonomy_df[level] != '') & (taxonomy_df[parent] != '')]
        n_relatives[level] = sub.groupby(parent)[level].nunique().to_dict()

    res = []
    for query, sub_df in tqdm(df.groupby('Query')):
        sub_df = sub_df.sort_values('Similarity', ascending=False)

        # Expected species = species of the query sequence itself
        expected_species = accession_species.get(query, '')
        if expected_species == '':
            continue
        expected = species_taxonomy.loc[expected_species]

        # Only queries whose own species is among the 100% hits
        top_hits = sub_df[sub_df['Similarity'] == 100]
        if expected_species not in top_hits['species'].values:
            continue

        # Query + its full taxonomy
        row = [query] + expected[tax_levels].tolist()
        for level in test_levels:
            # Skip levels without taxonomy (comparison would be meaningless)
            if expected[level] == '':
                row += [np.nan, False, np.nan]
                continue

            # Leave-one-out: best hit that does NOT belong to the expected taxon at this level
            other_hits = sub_df[sub_df[level] != expected[level]]
            if len(other_hits) > 0:
                similarity = other_hits['Similarity'].iloc[0]
            else:
                # All hits belong to the expected taxon: the closest other taxon is even further away
                similarity = sub_df['Similarity'].iloc[-1]

            # Isolated lineage: parent taxon has too few members (e.g. only one species in the genus)
            relatives = n_relatives[level].get(expected[parent_level[level]], 0)

            # Same taxon: closest hit of the expected taxon, but not the query sequence itself (species),
            # or of another lower taxon (e.g. genus: closest OTHER species of the same genus)
            same_hits = sub_df[(sub_df[level] == expected[level]) & (sub_df['Accession'] != query)]
            child = child_level[level]
            if child is not None:
                same_hits = same_hits[(same_hits[child] != expected[child]) & (same_hits[child] != '')]
            if len(same_hits) > 0:
                same_similarity = same_hits['Similarity'].iloc[0]
            else:
                same_similarity = np.nan

            row += [similarity, relatives >= min_relatives, same_similarity]
        res.append(row)

    # Taxonomy columns are capitalised (e.g. "Class"), similarity columns are lower case (e.g. "class")
    columns = ['Query'] + [i.capitalize() for i in tax_levels]
    for level in test_levels:
        columns += [level, f'{level} (has relatives)', f'{level} (same taxon)']
    res_df = pd.DataFrame(res, columns=columns)

    # 1) Thresholds for the whole database: optimal separation of "same taxon" and "other taxon" similarities
    summary, curves = summarise_levels(res_df, test_levels, fp_weight)
    if summary['n same taxon'].sum() == 0:
        print('Note: no "same taxon" hits found (e.g. only one sequence per species) - using the 95th percentile instead.')
    all_values = summary['Optimal threshold'].fillna(summary['95% other taxon'])
    thresholds_table = [['all', 'All', len(res_df)] + [round(i, 1) for i in all_values.values.tolist()]]
    print('Thresholds (All): ' + ','.join([str(i) for i in thresholds_table[0][3:]]))
    for _, row in summary.iterrows():
        print(f"  {row['Level']}: {row['Optimal threshold']} (false negatives {row['False negatives (%)']}%, "
              f"false positives {row['False positives (%)']}%)")

    # 2) Thresholds per taxon of the chosen level (e.g. per class)
    group_col = taxon_level.capitalize()
    group_summaries = []
    for taxon, group_df in res_df.groupby(group_col):
        if taxon == '' or len(group_df) < min_group_size:
            continue  # too few queries for a reliable threshold
        group_summary, _ = summarise_levels(group_df, test_levels, fp_weight, fallback=False)
        # Levels without enough data in this taxon: use the value of the whole database
        values = group_summary['Optimal threshold'].fillna(all_values)
        thresholds_table.append([taxon_level, taxon, len(group_df)] + [round(i, 1) for i in values.values.tolist()])
        group_summary.insert(0, group_col, taxon)
        group_summaries.append(group_summary)

    # One row per taxon ("All" first), one column per taxonomic level
    thresholds_df = pd.DataFrame(thresholds_table, columns=['Level', 'Taxon', 'n queries'] + test_levels)
    print(f'Thresholds calculated for {len(group_summaries)} taxa ({taxon_level}, min. {min_group_size} queries).')

    with pd.ExcelWriter(output_xlsx, engine='xlsxwriter') as writer:
        thresholds_df.to_excel(writer, sheet_name='Thresholds', index=False)
        summary.to_excel(writer, sheet_name='Summary', index=False)
        if len(curves) > 0:
            error_curves = pd.concat(curves)
            error_curves.to_excel(writer, sheet_name='Error curves', index=False)
        if len(group_summaries) > 0:
            pd.concat(group_summaries).to_excel(writer, sheet_name=f'Summary per {taxon_level}', index=False)
        res_df.to_excel(writer, sheet_name='Results', index=False)

    print(f'Plotting error curves...')
    plot_barcode_gap(error_curves, thresholds_df, summary, output_xlsx)

    return summary

def plot_barcode_gap(error_curves, thresholds_df, summary, output_xlsx):
    levels = error_curves['Level'].unique()
    for i, level in enumerate(levels):
        output_pdf = output_xlsx.parent / f'{output_xlsx.stem}_{i}_{level}.pdf'
        title = level.capitalize()
        fig = go.Figure()
        sub_df = error_curves[error_curves['Level'] == level]
        # Barcode Gap
        x_values = sub_df['Threshold'].values.tolist()
        y1_values = sub_df['False positives (%)'].values.tolist()
        y2_values = sub_df['False negatives (%)'].values.tolist()
        y3_values = sub_df['Weighted error'].values.tolist()
        fig.add_trace(go.Scatter(x=x_values, y=y1_values, name='False positives (%)', fill='tozeroy', mode='lines', line=dict(width=1), marker_color='Yellow'))
        fig.add_trace(go.Scatter(x=x_values, y=y2_values, name='False negatives (%)', fill='tozeroy', mode='lines', line=dict(width=1), marker_color='Red'))
        fig.add_trace(go.Scatter(x=x_values, y=y3_values, name='Weighted error', mode='lines', line=dict(width=1, dash='dot'), marker_color='Black'))
        # Threshold
        sub_df.sort_values(by=['Weighted error'], inplace=True)
        x_pos = summary.loc[summary['Level'] == level, 'Optimal threshold'].values[0]
        y_pos = sub_df.loc[sub_df['Threshold'] == x_pos, 'Weighted error'].values[0]
        fig.add_annotation(x=x_pos, y=y_pos, text=f'{x_pos:.1f}%', showarrow=True, ax=20, ay=20, xanchor='left', yanchor='middle', font=dict(size=10), arrowhead=2, arrowsize=1, arrowwidth=1)
        # Layout & Export
        fig.update_layout(template='simple_white', width=500, height=400, title=title)
        fig.update_layout(legend=dict( orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1))
        fig.update_xaxes(title_text='Similarity threshold (%)', autorange=False, showgrid=True, dtick=2, range=(100, 50))
        fig.update_yaxes(title_text='Share (%)', showgrid=True, dtick=20, range=(0,100))
        fig.write_image(output_pdf)


def add_errors(seq, sub_rate=0.01, ins_rate=0.0, del_rate=0.0, seed=None):
    random.seed(seed)
    bases = 'ACGT'
    new_seq = ''
    for base in seq:
        r = random.random()
        if r < del_rate:
            continue                                   # deletion: skip this base
        if r < del_rate + sub_rate:
            base = random.choice([b for b in bases if b != base])  # substitution
        new_seq += base
        if random.random() < ins_rate:
            new_seq += random.choice(bases)            # insertion after this base
    return new_seq

def create_test_fastas(work_folder, error_values, n_test=1000):
    # Random subset of the query sequences, one fasta per error rate
    records = list(SeqIO.parse(work_folder / 'tmp_query.fasta', 'fasta'))
    random.seed(1)
    records = random.sample(records, min(n_test, len(records)))

    test_fastas = {}
    for error in error_values:
        test_fasta = work_folder / f'tmp_test_{error}.fasta'
        f = open(test_fasta, 'w')
        for i, record in enumerate(records):
            seq = add_errors(str(record.seq), sub_rate=error, seed=i)  # seed per sequence = reproducible
            f.write(f'>{record.id}\n{seq}\n')
        f.close()
        test_fastas[error] = test_fasta
    return test_fastas

def evaluate_run(taxonomy_table, accession_species):
    # 1 = assigned to the expected species, 0 = not (other species, higher rank only, or no match)
    if not taxonomy_table.exists():
        return np.nan, np.nan
    taxonomy_df = pd.read_excel(taxonomy_table).fillna('')
    evaluation = []
    for _, row in taxonomy_df.iterrows():
        expected_species = accession_species.get(row['unique ID'], '')
        if expected_species == '':
            continue
        if row['Species'] == expected_species:
            evaluation.append(1)
        else:
            evaluation.append(0)
    if len(evaluation) == 0:
        return np.nan, np.nan
    return np.mean(evaluation), np.std(evaluation)

def dbtester(db_path, work_folder, n_cores, n_test):
    # Settings
    tasks = ['blastn', 'dc-megablast', 'megablast']
    max_target_seqs = [100]
    filter_modes = [1, 2]
    sim_ranges = [0, 0.5, 1]
    # Thresholds to compare: "All" of db_thresholds.xlsx, the default list, and the taxon-specific thresholds
    threshold_df = pd.read_excel(db_path / 'db_thresholds.xlsx', sheet_name='Thresholds')
    all_row = threshold_df[threshold_df['Taxon'] == 'All'].iloc[0]
    all_thresholds = ','.join([str(all_row[level]) for level in ['species', 'genus', 'family', 'order', 'class']])
    thresholds = [f'-thresholds {all_thresholds}', '-thresholds 98,95,90,87,85']
    taxon_levels = [i for i in threshold_df['Level'].unique() if i != 'all']
    if len(taxon_levels) > 0:
        thresholds.append(f'-threshold_level {taxon_levels[0]}')
    error_values = [0.0, 0.01, 0.02, 0.03]  # sub_rate

    # Test sequences (one fasta per error rate)
    test_fastas = create_test_fastas(work_folder, error_values, n_test=n_test)

    # Species lookup
    taxonomy_df = pd.read_parquet(db_path / 'db_taxonomy.parquet.snappy').fillna('')
    accession_species = dict(zip(taxonomy_df['Accession'], taxonomy_df['species']))

    res = []
    settings = list(product(error_values, tasks, max_target_seqs, filter_modes, sim_ranges, thresholds))
    print(f'{len(settings)} runs')
    for i, (error, task, max_target_seq, filter_mode, sim_range, threshold) in enumerate(tqdm(settings)):
        out = work_folder / f'run_{i}'
        taxonomy_table = out / f'run_{i}_taxonomy.xlsx'

        # Only run if not done yet (lets you restart after a crash); finished runs are still evaluated
        if not taxonomy_table.exists():
            if out.exists():
                shutil.rmtree(out)  # incomplete run from a crash
            cmd = (f'apscale_blast -db "{db_path}" -q "{test_fastas[error]}" -o "{out}" -n_cores {n_cores} -task {task} '
                   f'-max_target_seqs {max_target_seq} -filter {filter_mode} -sim_range {sim_range} {threshold}')
            try:
                tqdm.write(cmd)
                subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            except subprocess.CalledProcessError as e:
                tqdm.write(f'run_{i} failed: {e.stderr.strip()[-300:]}')

        # Evaluate and remove the large intermediate files of this run (only the taxonomy table is kept)
        eval_mean, eval_std = evaluate_run(taxonomy_table, accession_species)
        if (out / 'subsets').exists():
            shutil.rmtree(out / 'subsets')
        res.append([i, error, task, max_target_seq, filter_mode, sim_range, threshold, eval_mean, eval_std])

    # Table of all runs and their settings (best runs first)
    res_df = pd.DataFrame(res, columns=['Run', 'Error_rate', 'Task', 'Max_target_seqs', 'Filter',
                                        'Sim_range', 'Threshold', 'Evaluation mean', 'Evaluation std'])
    res_df.sort_values('Evaluation mean', ascending=False, inplace=True)
    res_df.to_excel(work_folder / 'runs.xlsx', index=False)

def cleanup(work_folder):
    # Remove all temporary files (tmp_*); results (leave_one_out.xlsx, runs.xlsx, run_* folders) are kept
    for file in work_folder.glob('tmp_*'):
        file.unlink()

def main(db_path, n_cores, fullrun, n_test, taxon_level='class', fp_weight=1):
    db_path = Path(db_path)

    # All files are written here - the database folder itself is never changed
    work_folder = db_path.parent / 'dbbenchmark' / db_path.name
    os.makedirs(work_folder, exist_ok=True)

    # 1) Calculate thresholds
    blast_csv = blastn(db_path, work_folder, n_cores=n_cores)
    if blast_csv is None:
        return
    leave_one_out(blast_csv, db_path / 'db_taxonomy.parquet.snappy', db_path / 'db_thresholds.xlsx', taxon_level=taxon_level, fp_weight=fp_weight)

    # 2) Run all combinations of all settings
    if fullrun == True:
        dbtester(db_path, work_folder, n_cores, n_test)

    # 3) Remove temporary files
    cleanup(work_folder)

if __name__ == '__main__':
    # Direct call: python e_dboptimiser.py <database> <n_cores>
    main(sys.argv[1], int(sys.argv[2]), n_test=1000, fullrun=True)