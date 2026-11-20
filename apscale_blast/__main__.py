import argparse
import multiprocessing
import os
import sys
import datetime
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import subprocess
import importlib
import io
from ete3 import NCBITaxa, Tree
import sqlite3
from pandas import read_sql_query, read_sql_table
from tqdm import tqdm
import math
# APSCALE-blast
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter
from apscale_blast.d_gbif_validation import main as d_gbif_validation
from apscale_blast.e_dboptimiser import main as e_dboptimiser
# BOLDigger3
from boldigger3 import id_engine
from boldigger3 import add_metadata
from boldigger3 import select_top_hit
from boldigger3 import download_database

import warnings

warnings.filterwarnings("ignore", message="taxid .* was translated into")

ncbi = NCBITaxa()


def organism_filter(name):
    try:
        name = int(name)
    except ValueError:
        pass
    if isinstance(name, int):
        try:
            # Initialize NCBITaxa database
            ncbi = NCBITaxa()
            res = ncbi.get_taxid_translator([name])
            organism = f'{res[name]} (taxid:{name})'
            return organism
        except:
            return f'{name} not found!'
    else:
        try:
            # Initialize NCBITaxa database
            ncbi = NCBITaxa()
            res = ncbi.get_name_translator([name])
            taxid = res[name][-1]
            organism = f'{name} (taxid:{taxid})'
            return organism
        except:
            return f'{name} not found!'


def issue_list_flag(merged_df, issue_list):
    # Merge on shared columns
    merged_with_issues = merged_df.merge(
        issue_list,
        on=["Species", "Family", "Phylum"],
        how="left",  # keep all rows from merged_df
        suffixes=("", "_issue")  # avoid name conflicts
    )

    merged_with_issues = merged_with_issues.fillna('')
    col = "Status"  # column to move
    merged_with_issues = merged_with_issues[[c for c in merged_with_issues.columns if c != col] + [col]]

    return merged_with_issues


def merge_reblast(taxonomy_df, taxonomy_reblast_df, db_name):
    # Example key column (adjust as needed)
    key_col = "unique ID"

    # Merge both tables
    merged = taxonomy_df.merge(
        taxonomy_reblast_df,
        on=key_col,
        suffixes=("", "_reblast"),
        how="outer"  # keep all sequences
    )

    # Compare similarity columns
    def choose_row(row):
        # Case 1: Only one of the two has data
        if pd.isna(row["Similarity"]):
            row["Similarity"] = row["Similarity_reblast"]
            row["Reblast result"] = True
        elif pd.isna(row["Similarity_reblast"]):
            row["Reblast result"] = ''
        # Case 2: Both exist — compare
        elif row["Similarity_reblast"] > row["Similarity"]:
            for col in taxonomy_reblast_df.columns:
                if col != key_col:
                    row[col] = row[f"{col}_reblast"]
            row["Reblast result"] = True
        else:
            row["Reblast result"] = ''
        return row

    merged = merged.apply(choose_row, axis=1)

    # Drop the reblast columns (clean up)
    merged = merged[[c for c in merged.columns if not c.endswith("_reblast")]]

    return merged


def run_boldigger3(query_fasta):
    # Collect input values
    thresholds = [97, 95, 90, 85, 75, 50]
    query_fasta = str(query_fasta)
    mode = 3  # exhaustive Search
    db = 2  # ANIMAL SPECIES-LEVEL LIBRARY (PUBLIC + PRIVATE)

    # download the current metadata from BOLD
    metadata_download.main()

    # run the id engine
    id_engine.main(
        query_fasta,
        db,
        mode,
    )

    # add additional data via the metadata
    add_metadata.main(query_fasta)

    # select the top hit
    select_top_hit.main(query_fasta, thresholds)


def merge_boldigger3(taxonomy_df, boldigger_df, db_name):
    # Example key column (adjust as needed)
    key_col = "unique ID"

    boldigger_df.columns = [i.capitalize() for i in boldigger_df.columns.tolist()]
    boldigger_df = boldigger_df.rename(columns={"Id": "unique ID", "Pct_identity": "Similarity"})

    # Merge both tables
    merged = taxonomy_df.merge(
        boldigger_df,
        on=key_col,
        suffixes=("", "_reblast"),
        how="outer"  # keep all sequences
    )

    # Compare similarity columns
    def choose_row(row):
        # Case 1: Only one of the two has data
        if pd.isna(row["Similarity"]):
            row["Similarity"] = row["Similarity_reblast"]
            row["Reblast result"] = True
        elif pd.isna(row["Similarity_reblast"]):
            row["Reblast result"] = ''
        # Case 2: Both exist — compare
        elif row["Similarity_reblast"] > row["Similarity"]:
            for col in boldigger_df.columns:
                reblast_col = f"{col}_reblast"
                if col != key_col and reblast_col in merged.columns:
                    row[col] = row[reblast_col]
            row["Reblast result"] = True
        else:
            row["Reblast result"] = ''
        return row

    merged = merged.apply(choose_row, axis=1)
    merged = merged.fillna('')

    # Drop the reblast columns (clean up)
    merged = merged[[c for c in merged.columns if not c.endswith("_reblast")]]

    return merged


def check_blast_db(db_path):
    print('{}: Validating database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
    try:
        subprocess.run(
            ["blastdbcmd", "-db", str(db_path), "-info"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True
        )
        print('{}: Your database looks fine.'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        return True
    except:
        print('{}: Database could not be verified!'.format(datetime.datetime.now().strftime('%H:%M:%S')))
        return False


def convert_ncbi_db(db_path):
    print('{}: Converting ncbi database.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    # 1) accession + taxid for every sequence in the db
    out = subprocess.run(["blastdbcmd", "-db", str(db_path / db_path.name), "-entry", "all", "-outfmt", "%a\t%T"],
                         capture_output=True, text=True, check=True).stdout
    acc = pd.read_csv(io.StringIO(out), sep="\t", names=["Accession", "taxid"])

    # 2) taxid -> ranks
    ncbi = NCBITaxa()
    RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

    def taxids_to_df(taxids, ranks):
        taxids = pd.Series(taxids).dropna().astype(int).unique()
        lineages = {}
        for t in tqdm(taxids, desc="Converting taxids to lineages"):
            try:
                lineages[t] = ncbi.get_lineage(t)
            except ValueError:
                lineages[t] = []
        nodes = {n for lin in lineages.values() for n in lin}
        rank_of = ncbi.get_rank(nodes)
        name_of = ncbi.get_taxid_translator(nodes)
        rows = []
        for t, lin in lineages.items():
            r2n = {rank_of[n]: name_of[n] for n in lin}
            rows.append([t] + [r2n.get(r) for r in ranks])
        return pd.DataFrame(rows, columns=["taxid"] + ranks)

    tax = taxids_to_df(acc["taxid"], RANKS)

    # 3) final table: Accession | Superkingdom | ... | Species
    taxonomy_df = (acc.merge(tax, on="taxid", how="left").drop(columns="taxid").rename(columns=str.capitalize))
    taxonomy_df.fillna('', inplace=True)
    taxonomy_df.rename(columns={"Kingdom": "superkingdom"}, inplace=True)
    taxonomy_df.columns = [i.lower() if i != 'Accession' else i.capitalize() for i in taxonomy_df.columns]
    taxonomy_file = db_path / "db_taxonomy.parquet.snappy"
    taxonomy_df.to_parquet(taxonomy_file, index=False)

    files = ['taxdb.btd', 'taxdb.bti', 'taxonomy4blast.sqlite3']
    for file in files:
        file_path = db_path / file
        if file_path.exists():
            os.unlink(file_path)

    # rebuild db under the new name (renaming v5 files breaks the LMDB index)
    old_db = str(db_path / db_path.name)
    new_db = str(db_path / "db")
    fasta = db_path / "tmp_seqs.fasta"

    with open(fasta, "w") as fh:
        subprocess.run(["blastdbcmd", "-db", old_db, "-entry", "all", "-outfmt", "%f"], stdout=fh, check=True)
    subprocess.run(["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", new_db, "-title", "db"], check=True)

    # remove old db files and temp files
    for f in db_path.glob(f"{db_path.name}.*"):
        f.unlink()
    fasta.unlink()


def get_package_versions(pkg):
    try:
        version = importlib.metadata.version(pkg)
        return version
    except importlib.metadata.PackageNotFoundError:
        return


def settings_to_df(args):
    settings = {
        # Main settings
        'database': args.database,
        'query_fasta': args.query_fasta,
        'out': args.out,
        # BLASTn settings
        'n_cores': args.n_cores,
        'task': args.task,
        'subset_size': args.subset_size,
        'max_target_seqs': args.max_target_seqs,
        # Filter settings
        'thresholds': args.thresholds,
        'threshold_level': args.threshold_level,
        'filter': args.filter,
        'masking': args.masking,
        'rating_range': args.rating_range,
        'sim_range': args.sim_range,
        # Re-BLAST settings
        'reblast_db': args.reblast_db,
        'reblast_sim': args.reblast_sim,
        'reblast_only': args.reblast_only,
        # Misc
        'blastn_exe': args.blastn_exe,
    }

    df = pd.DataFrame(settings.items(), columns=['Setting', 'Value'])
    return df


def main():
    """
    APSCALE-blastn suite
    Command-line tool to run and filter BLASTn searches.
    """

    apscale_blast_version = get_package_versions('apscale_blast')

    # Introductory message with usage examples
    message = f"""
    APSCALE-blast command line tool - {apscale_blast_version}

    Example commands:
      $ apscale_blast -h

    Basic blast 
      $ apscale_blast -db ./MIDORI2_UNIQ_NUC_GB266_srRNA_BLAST -q ./12S_apscale_ESVs.fasta -o ./12S_apscale_ESV

    GBIF occurrence module (latitude = 50.6, longitude = 6.4, radius = 200km) 
      $ apscale_blast gbif -latlonrad 50.6,6.4,200 -q ./12S_apscale_ESV/12S_apscale_ESV_blastn.xlsx

    GBIF occurrence module (countries: Germany, Luxembourg, France) 
    $ apscale_blast gbif -countries DE,LU,FR -q ./12S_apscale_ESV/12S_apscale_ESV_blastn.xlsx

    GBIF synonym check only (no occurrence search)
      $ apscale_blast gbif -synonyms -q ./12S_apscale_ESV/12S_apscale_ESV_blastn.xlsx

    Re-BLAST of low-similarity sequences against a second database
      $ apscale_blast -db ./MIDORI2_UNIQ_NUC_GB266_srRNA_BLAST -db2 ./my_local_12S_db -q ./12S_apscale_ESVs.fasta -o ./12S_apscale_ESV

    Re-BLAST of an existing run (no new BLAST against -db)
      $ apscale_blast -reblast_only ./12S_apscale_ESV/12S_apscale_ESV_taxonomy.xlsx -db2 ./my_local_12S_db -q ./12S_apscale_ESVs.fasta

    Database creator (streamlit GUI)
      $ apscale_blast dbcreator

    Database optimiser (leave-one-out thresholds + settings benchmark)
      $ apscale_blast dboptimiser -db ./MIDORI2_UNIQ_NUC_GB266_srRNA_BLAST -n_cores 8

    """
    print(message)

    ####################################################################################################################
    ## 1) dbcreator
    if len(sys.argv) > 1 and sys.argv[1] == "dbcreator":
        # Determine theme
        theme = "dark"  # default
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
            script_path = os.path.join(script_dir, 'c_dbcreator.py')

            print(f"Starting APSCALE-blast dbcreator in {theme} mode. Press CTRL + C to exit.")
            subprocess.run([
                'streamlit', 'run', script_path,
                '--theme.base', theme,
                '--server.address', 'localhost'
            ])
        except KeyboardInterrupt:
            print("Exiting...")
            sys.exit()

    ####################################################################################################################
    ## 2) gbif module
    if len(sys.argv) > 1 and sys.argv[1] == "gbif":
        print("Starting APSCALE-blast GBIF module.")
        print('Each species will be validated against the GBIF occurrence database.')

        parser = argparse.ArgumentParser(
            prog="APSCALE-blast GBIF module",
            description=f'APSCALE-blast {apscale_blast_version}'
        )

        main_settings = parser.add_argument_group("Main settings")
        main_settings.add_argument('-latlonrad', type=str, required=False, default=False,
                                   help='Latitude, longitude, and radius to build the search polygon (example: "50.1,6.2,100")')
        main_settings.add_argument('-continent', type=str, required=False, default=False,
                                   help='Continent to search (e.g., "europe"')
        main_settings.add_argument('-countries', type=str, required=False, default=False,
                                   help='Countries to search (e.g., "DE,LU,FR"')
        main_settings.add_argument('-similarity', '-s', type=float, required=False, default=98,
                                   help='Minimum similarity score to check for occurrences.')
        main_settings.add_argument('-habitat', action='store_true',
                                   help='Fetch the habitat of each species from GBIF.')
        main_settings.add_argument('-query', '-q', type=str, required=True,
                                   help='APSCALE-blast taxonomy table.')
        main_settings.add_argument('-f7', '-flag7', action='store_true',
                                   help='Enable flag 7 (species occurrences of the same genus).')
        main_settings.add_argument('-synonyms', action='store_true',
                                   help='Only replace synonyms by their accepted names (no occurrence search, no search area needed).')
        args = parser.parse_args(sys.argv[2:])

        if args.f7 == False and args.synonyms == False:
            print(' Note: Flag 7 is disabled (species occurrences of the same genus).')
        if args.habitat == False:
            print(' Note: Habitat search is disabled.')

        d_gbif_validation(Path(args.query), args.continent, args.countries, args.latlonrad, args.similarity,
                          args.habitat, args.f7, synonyms_only=args.synonyms)
        print(f'\n{datetime.datetime.now().strftime("%H:%M:%S")}: Done! Have a nice day!')
        return

    ####################################################################################################################
    ## 3) dboptimiser module
    if len(sys.argv) > 1 and sys.argv[1] == "dboptimiser":
        print("Starting APSCALE-blast database optimiser.")
        print('Similarity thresholds are estimated from the database using a leave-one-out approach.')

        parser = argparse.ArgumentParser(
            prog="APSCALE-blast dboptimiser module",
            description=f'APSCALE-blast {apscale_blast_version}'
        )

        main_settings = parser.add_argument_group("Main settings")
        main_settings.add_argument('-database', '-db', type=str, required=True,
                                   help='Path to the local reference database.')
        main_settings.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 2,
                                   help='Number of CPU cores to use. [Default: CPU count - 2]')
        main_settings.add_argument('-n_test', type=int, default=1000,
                                   help='Number of test sequences for the settings benchmark. [Default: 1000]')
        main_settings.add_argument('-full_benchmark', '-f', action='store_true',
                                   help='Calculate the thresholds and the settings benchmark.')
        main_settings.add_argument('-taxon_level', type=str, default='class',
                                   choices=['phylum', 'class', 'order', 'family', 'genus'],
                                   help='Taxonomic level for taxon-specific thresholds. [Default: class]')
        main_settings.add_argument('-fp_weight', type=float, default=1,
                                   help='Weight of false positives (wrong assignments) relative to false negatives '
                                        '(too high-level assignments). Values > 1 give more conservative thresholds. [Default: 1]')
        args = parser.parse_args(sys.argv[2:])

        # Verify the database before starting
        database = Path(args.database.strip('"'))
        if check_blast_db(database / 'db') == False:
            return

        if args.full_benchmark == True:
            full_benchmark = True
        else:
            full_benchmark = False
            print(' Note: Settings benchmark is disabled (thresholds only).')

        e_dboptimiser(database, args.n_cores, full_benchmark, args.n_test, args.taxon_level, args.fp_weight)
        print(f'\n{datetime.datetime.now().strftime("%H:%M:%S")}: Done! Have a nice day!')
        return

    ####################################################################################################################
    ## 4) APSCALE-blast
    parser = argparse.ArgumentParser(description=f'APSCALE-blast {apscale_blast_version}')

    # === Main settings ===
    main_settings = parser.add_argument_group("Main settings")
    main_settings.add_argument('-database', '-db', type=str, required=False,
                               help='Path to the local reference database (not needed with -reblast_only).')
    main_settings.add_argument('-query_fasta', '-q', type=str, required=True,
                               help='Path to the query FASTA file.')
    main_settings.add_argument('-out', '-o', type=str, default='./blastn',
                               help='Output directory. A new folder will be created here. [Default: ./blastn]')

    # === BLASTn settings ===
    blastn_settings = parser.add_argument_group("BLASTn settings")
    blastn_settings.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 2,
                                 help='Number of CPU cores to use. [Default: CPU count - 2]')
    blastn_settings.add_argument('-task', type=str, default='megablast',
                                 help='BLASTn task: blastn, megablast, or dc-megablast. [Default: megablast]')
    blastn_settings.add_argument('-subset_size', type=int, default=100,
                                 help='Number of sequences per FASTA subset. [Default: 100]')
    blastn_settings.add_argument('-max_target_seqs', type=int, default=20,
                                 help='Maximum number of hits retained per query. [Default: 20]')

    # === Filter settings ===
    filter_settings = parser.add_argument_group("Filter settings")
    filter_settings.add_argument('-thresholds', type=str, default=None,
                                 help='Similarity thresholds for species, genus, family, order, class (comma-separated, e.g. "98,95,90,87,85"). '
                                      'If set, these are used instead of db_thresholds.xlsx. '
                                      '[Default: db_thresholds.xlsx of the database, otherwise 98,95,90,87,85]')
    filter_settings.add_argument('-threshold_level', type=str, default='all',
                                 choices=['all', 'phylum', 'class', 'order', 'family', 'genus'],
                                 help='Thresholds of db_thresholds.xlsx: "all", or taxon-specific (e.g. "class"). '
                                      'Taxa based on < 50 sequences use "all". [Default: all]')
    filter_settings.add_argument('-filter', type=str, default='2',
                                 help='Filtering mode: 1 = e-value → similarity, 2 = similarity → e-value, 3 = similarity [Default: 2]')
    filter_settings.add_argument('-masking', action='store_false',
                                 help='Disable masking. Masking is enabled by default.')
    filter_settings.add_argument('-rating_range', type=int, default=5,
                                 help='Range of allowed rating values. [Default: 5]')
    filter_settings.add_argument('-sim_range', type=int, default=0,
                                 help='Value that is subtracted from the highest similarity to retain hits. [Default: 0]')

    # === Re-BLAST settings ===
    reblast_settings = parser.add_argument_group("Re-BLAST settings")
    reblast_settings.add_argument('-reblast_db', '-db2', type=str,
                                  help='Path to an optional secondary database for re-BLAST. Choose "boldigger" for boldigger3 assignment.')
    reblast_settings.add_argument('-reblast_sim', type=float, default=99,
                                  help='Sequences below this similarity are re-BLASTed against db2. [Default: 99]')
    reblast_settings.add_argument('-reblast_only', type=str, default=None,
                                  help='Path to an existing APSCALE-blast taxonomy table (.xlsx) that is re-BLASTed against -db2. '
                                       'The initial BLAST is skipped and -db is not needed.')

    # === Misc settings ===
    misc_settings = parser.add_argument_group("Misc")
    misc_settings.add_argument('-blastn_exe', type=str, default='blastn',
                               help='Path to the BLASTn executable. [Default: blastn]')

    # Parse the arguments
    args = parser.parse_args()

    # Handle missing arguments interactively for both commands
    if not args.database and not args.query_fasta:
        args.database = input("Please enter PATH to database: ").strip('"')
        args.query_fasta = input("Please enter PATH to query fasta: ").strip('"')

        # Set output directory if default value is used
        if args.out == './':
            args.out = str(args.query_fasta).replace('.fasta', '')
            if not os.path.isdir(args.out):
                os.mkdir(Path(args.out))  # Create the output directory

    ## CHECK IF FILES ALREADY EXIST
    project_folder = Path(args.out)  # Use the output directory specified by the user

    # Check the required inputs
    if args.reblast_only:
        # Re-BLAST of an existing taxonomy table: only -q, -db2 and -reblast_only are needed
        reblast_only_xlsx = Path(args.reblast_only.strip('"'))
        if not reblast_only_xlsx.is_file():
            print(f'\nError: Could not find the taxonomy table: {reblast_only_xlsx}')
            return
        if not args.reblast_db:
            print('\nError: Please provide a secondary database (-db2) for the re-BLAST.')
            return
    elif not args.database:
        print('\nError: Please provide a database (-db).')
        return

    # Filter mode (used for the BLAST and the re-BLAST)
    if not args.filter:
        filter_mode = "2"
    else:
        filter_mode = str(args.filter)

    if not args.reblast_only:
        # Convert db to Path
        database = Path(args.database.strip('"'))
        db_name = Path(database).stem

        # Verify db
        if ' ' in str(database):
            print('\nError: Your database PATH contains white spaces! Please move your database to a different folder!')
            print('')
            return
        db_path = database / 'db'
        # First check if a blast+ GenBank db can be converted
        taxonomy4blast_sqlite3 = db_path.parent / 'taxonomy4blast.sqlite3'
        if taxonomy4blast_sqlite3.exists():
            ncbi_db_path = db_path.parent
            convert_ncbi_db(ncbi_db_path)
        # Then verify the database
        db_check = check_blast_db(db_path)
        if db_check == False:
            print('')
            return

        # Run the BLASTn function
        a_blastn(args.blastn_exe,
                 args.query_fasta.strip('"'),
                 database,
                 project_folder,
                 args.n_cores,
                 args.task,
                 args.subset_size,
                 args.max_target_seqs,
                 args.masking,
                 )

        # Handle the 'filter' command
        if not os.path.isfile(Path(project_folder).joinpath('log.txt')):
            print('\nError: Could not find the BLAST results folder!')
            return
        print('')
        if filter_mode == "1":
            print(
                f'{datetime.datetime.now().strftime("%H:%M:%S")}: Using filter mode 1: similarity (-{args.sim_range}%) -> e-value -> rating (optional)')
        elif filter_mode == "2":
            print(
                f'{datetime.datetime.now().strftime("%H:%M:%S")}: Using filter mode 2: e-value -> similarity (-{args.sim_range}%) -> rating (optional)')
        else:
            print(
                f'{datetime.datetime.now().strftime("%H:%M:%S")}: Using filter mode 3: similarity (-{args.sim_range}%) -> rating (optional)')
        b_filter(project_folder, database, args.thresholds, args.n_cores, filter_mode, args.rating_range,
                 args.sim_range, args.threshold_level)

        # Taxonomy table for the re-BLAST (if defined) and the merged output
        name = project_folder.name
        taxonomy_xlsx = Path('{}/{}_taxonomy.xlsx'.format(project_folder, name))
        reblast_output_xlsx = Path('{}/{}_reblast_taxonomy.xlsx'.format(project_folder, name))

    else:
        # Existing taxonomy table: all re-BLAST results are written next to it
        taxonomy_xlsx = reblast_only_xlsx
        project_folder = taxonomy_xlsx.parent
        db_name = taxonomy_xlsx.stem
        reblast_output_xlsx = taxonomy_xlsx.parent / taxonomy_xlsx.name.replace('.xlsx', '_reblast.xlsx')
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Re-BLAST of the existing taxonomy table "{taxonomy_xlsx.name}".')

    # Perform reblast if defined
    if args.reblast_db:
        print(f'\n{datetime.datetime.now().strftime("%H:%M:%S")}: Performing a re-blast of low-smiliarity sequences.')

        # Define required folders and files (one folder per secondary database, so different -db2 runs don't mix)
        reblast_db = Path(args.reblast_db.strip('"'))
        reblast_folder = project_folder / f'reblast_{reblast_db.name}'
        os.makedirs(reblast_folder, exist_ok=True)
        query_fasta = Path(args.query_fasta.strip('"'))
        low_similarity_query_fasta = reblast_folder / 'lowsim.fasta'

        # Extract all sequences without a proper assignment and write them to a new fasta:
        # below the similarity threshold (default: 99%), without a species-level assignment,
        # or ambiguous within a genus (F2: "Genus a/b", F3: "Genus sp.")
        # (sequences without any BLAST hit have a similarity of 0 and are included)
        reblast_sim = args.reblast_sim
        taxonomy_df = pd.read_excel(taxonomy_xlsx).fillna('')
        low_similarity = taxonomy_df['Similarity'] < reblast_sim
        no_species = taxonomy_df['Species'] == ''
        flag_f2_f3 = (taxonomy_df['Flag'].astype(str).str.contains('F2 (', regex=False) |
                      taxonomy_df['Flag'].astype(str).str.contains('F3 (', regex=False))
        lowsim_df = taxonomy_df[low_similarity | no_species | flag_f2_f3]
        lowsim_hashes = set(lowsim_df['unique ID'].values.tolist())
        n_lowsim = len(lowsim_hashes)
        print(
            f'{datetime.datetime.now().strftime("%H:%M:%S")}: Found {n_lowsim} sequences for re-assignment (<{reblast_sim}% similarity, no species-level assignment, or F2/F3).')

        with open(low_similarity_query_fasta, "w") as out_handle:
            for record in SeqIO.parse(query_fasta, "fasta"):
                if record.id in lowsim_hashes:
                    # Write to output FASTA
                    SeqIO.write(record, out_handle, "fasta")

        if n_lowsim == 0:
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Nothing to re-BLAST.')

        elif reblast_db == Path('boldigger'):
            pass
            # run_boldigger3(low_similarity_query_fasta)
            # boldigger_results_xlsx = reblast_folder / 'boldigger3_data' / 'lowsim_identification_result.xlsx'
            # boldigger_df = pd.read_excel(boldigger_results_xlsx).fillna('')
            # merged_df = merge_boldigger3(taxonomy_df, boldigger_df, 'boldigger3')
            # col = "Status"  # column to move
            # merged_df = merged_df[[c for c in merged_df.columns if c != col] + [col]]
            # merged_df.to_excel(reblast_output_xlsx, index=False)
            # print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Finished re-assignment using BOLDigger3!')

        elif check_blast_db(reblast_db / 'db') == False:
            print('\nError: The secondary database (-db2) could not be verified!')

        else:
            # Run the BLASTn function
            a_blastn(args.blastn_exe,
                     low_similarity_query_fasta,
                     reblast_db,
                     reblast_folder,
                     args.n_cores,
                     args.task,
                     args.subset_size,
                     args.max_target_seqs,
                     args.masking,
                     )
            # Run the filter command
            print('')
            b_filter(reblast_folder, reblast_db, args.thresholds, args.n_cores, filter_mode, args.rating_range,
                     args.sim_range, args.threshold_level)

            taxonomy_reblast_xlsx = reblast_folder / f'{reblast_folder.name}_taxonomy.xlsx'
            taxonomy_reblast_df = pd.read_excel(taxonomy_reblast_xlsx).fillna('')
            merged_df = merge_reblast(taxonomy_df, taxonomy_reblast_df, db_name)
            col = "Status"  # column to move
            merged_df = merged_df[[c for c in merged_df.columns if c != col] + [col]]
            merged_df.to_excel(reblast_output_xlsx, index=False)
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Finished re-blast! Saved as "{reblast_output_xlsx.name}".')

    # If defined, use an issue list
    issue_list_xlsx = False
    if issue_list_xlsx == True:
        issue_list_flag(taxonomy_df, issue_list_df)

    # save settings to df
    df = settings_to_df(args)
    if args.reblast_only:
        settings_xlsx = Path('{}/{}_reblast_run_settings.xlsx'.format(project_folder, db_name))
    else:
        settings_xlsx = Path('{}/{}_blastn_run_settings.xlsx'.format(project_folder, db_name))
    df.to_excel(settings_xlsx, index=False)

    ## Finish script
    print(f'\n{datetime.datetime.now().strftime("%H:%M:%S")}: Done! Have a nice day!')


# Run the main function if script is called directly
if __name__ == "__main__":
    main()