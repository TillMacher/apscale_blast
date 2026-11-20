import argparse
import os
from pathlib import Path
from Bio import SeqIO
import datetime
import subprocess
import shutil
import time
import glob
import multiprocessing
import threading
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import pyarrow as pa
import pyarrow.parquet as pq
import sys
import os
from playwright.sync_api import sync_playwright
import time
from Bio import SeqIO
import json
import pandas as pd
from ete3 import NCBITaxa, Tree
import requests
import xmltodict
import random
import tkinter as tk
from tkinter import simpledialog, messagebox
from tqdm import tqdm
from itertools import islice
from Bio.SeqIO.FastaIO import SimpleFastaParser

# Lock for thread-safe print statements
print_lock = threading.Lock()

# Ensure multiprocessing works when frozen (e.g. for executables)
multiprocessing.freeze_support()

def fasta_subset(fasta_file, subset_size):
    """
    Splits a large fasta file into smaller subsets for parallel processing.
    """
    print(f'{datetime.datetime.now():%H:%M:%S}: Creating subset(s) from fasta file.')

    subset_size = int(subset_size)
    fasta_file = Path(fasta_file)

    # Create a fresh directory for the subsets
    subset_folder = fasta_file.parent / 'fasta_subsets'
    subset_folder.mkdir(parents=True, exist_ok=True)
    for f in subset_folder.glob('*.fasta'):
        f.unlink()

    # Read subset_size sequences at a time and write each chunk in one go
    n_subsets = 0
    with open(fasta_file) as handle:
        records = SimpleFastaParser(handle)
        while chunk := list(islice(records, subset_size)):
            n_subsets += 1
            with open(subset_folder / f'subset_{n_subsets}.fasta', 'w') as out:
                out.writelines(f'>{header}\n{seq}\n' for header, seq in chunk)

    print(f'{datetime.datetime.now():%H:%M:%S}: Created {n_subsets} subset(s) from fasta file.')
    return subset_folder

def json_to_csv(json_data, blastn_json_path):
    """
    Converts BLAST JSON output into a table with specific columns, including e-value and percentage identity.

    Args:
        json_data (dict): Parsed JSON data from BLAST output.

    Returns:
        pd.DataFrame: A DataFrame with the extracted information.
    """
    # Extract the relevant portion of the JSON data
    results = json_data["BlastOutput2"]
    table_data = []

    # Iterate over each report in the JSON
    for result in results:
        query_title = result["report"]["results"]["search"]["query_title"]
        hits = result["report"]["results"]["search"]["hits"]

        # Iterate over hits to collect data
        for hit in hits:
            try:
                for desc in hit["description"]:
                    taxid = desc["taxid"]

                    # Get the first HSP for e-value and percentage identity
                    if hit["hsps"]:
                        hsp = hit["hsps"][0]
                        evalue = hsp["evalue"]
                        identity = hsp["identity"]
                        align_len = hsp["align_len"]
                        percentage_identity = (identity / align_len) * 100
                    else:
                        evalue = None
                        percentage_identity = None

                    # Append to the table data
                    table_data.append({
                        "unique ID": query_title,
                        "Sequence ID": taxid,
                        "Similarity": percentage_identity,
                        "evalue": evalue,
                    })
            except KeyError:
                # Append to the table data
                table_data.append({
                    "unique ID": 'No Match',
                    "Sequence ID": 'No Match',
                    "Similarity": 0,
                    "evalue": 1,
                })

    # Convert to a pandas DataFrame
    df = pd.DataFrame(table_data)

    # Write to csv file
    blastn_csv_path = Path(str(blastn_json_path).replace('.json', '.csv'))

    # Cannot use two separators - ugly workaround required to match blast output format with two seps!
    # Save with a single-character delimiter
    temp_path = blastn_csv_path.with_suffix(".temp.csv")
    df.to_csv(temp_path, sep=";", index=False, header=False)
    # Replace the single-character delimiter with a multi-character one
    with open(temp_path, "r") as temp_file, open(blastn_csv_path, "w") as final_file:
        content = temp_file.read()
        final_file.write(content.replace(";", ";;"))
    # Optionally delete the temporary file
    os.remove(temp_path)
    os.remove(blastn_json_path)

def blastn_parallel(fasta_file, n_subsets, blastn_subset_folder, blastn_exe, db_folder, i, print_lock, task, max_target_seqs, masking):
    """
    Runs a single BLASTN job on a subset of the fasta file.

    Args:
        fasta_file (str): Path to the subset fasta file.
        n_subsets (int): Total number of subsets.
        blastn_subset_folder (Path): Folder to store BLASTN output.
        blastn_exe (str): Path to the BLASTN executable.
        db_folder (Path): Path to the BLASTN database.
        i (int): Subset index.
        print_lock (threading.Lock): Lock for synchronized printing.
        task (str): BLAST task (e.g., 'megablast', 'blastn').
        max_target_seqs (int): Maximum target sequences to report.
    """
    blastn_csv = blastn_subset_folder.joinpath(Path(fasta_file).stem + '_' + task + '.csv')

    # Skip if output already exists
    if os.path.isfile(blastn_csv) and os.path.getsize(blastn_csv) > 0:
        with print_lock:
            print('{}: Skipping {} (already exists and is not empty).'.format(datetime.datetime.now().strftime('%H:%M:%S'),
                                                             blastn_csv.stem))
        time.sleep(1)
    elif masking == False:
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Masking is disabled!')
        # Run the BLASTN command
        subprocess.call([blastn_exe, '-task', task, '-db', str(db_folder), '-query', str(fasta_file),
                         '-num_threads', str(1), '-max_target_seqs', str(max_target_seqs),
                         '-dust', 'no', '-soft_masking', 'false',
                         '-outfmt', '6 delim=;; qseqid sseqid pident evalue', '-out', str(blastn_csv)])
        with print_lock:
            print('{}: Finished tax. assignment for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1,
                                                                 n_subsets))

    else:
        # Run the BLASTN command
        subprocess.call([blastn_exe, '-task', task, '-db', str(db_folder), '-query', str(fasta_file),
                         '-num_threads', str(1), '-max_target_seqs', str(max_target_seqs),
                         '-outfmt', '6 delim=;; qseqid sseqid pident evalue tscore', '-out', str(blastn_csv)])
        with print_lock:
            print('{}: Finished tax. assignment for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1,
                                                                 n_subsets))

def main(blastn_exe, query_fasta, blastn_database, project_folder, n_cores, task, subset_size, max_target_seqs, masking):
    """
    Improved BLASTN function that utilizes multithreading for faster performance.

    Args:
        blastn_exe (str): Path to the BLASTN executable.
        query_fasta (str): Path to the input fasta file.
        blastn_database (str): Path to the BLASTN database.
        project_folder (str): Path to the project directory for saving output.
        n_cores (int): Number of cores to use for parallel processing.
        task (str): BLAST task (e.g., 'megablast', 'blastn').
        subset_size (int): Size of fasta file subsets.
        max_target_seqs (int): Maximum target sequences to report.
    """
    # Check if the accession table was already written
    filename = Path(query_fasta).stem.replace('.', '_').replace(' ', '_')
    accession_table = project_folder.joinpath('{}.parquet.snappy'.format(filename))

    if accession_table.exists():
        print(f'{datetime.datetime.now().strftime('%H:%M:%S')}: Skipping taxonomic assignment (already done).')
        return

    # Split fasta file into subsets
    subset_folder = fasta_subset(query_fasta, subset_size)

    project_folder = Path(project_folder)
    fasta_files = sorted(glob.glob(str(subset_folder) + '/*.fasta'))
    n_subsets = len(fasta_files)

    # Map task names to valid BLAST task identifiers
    task_mapping = {
        'Highly similar sequences (megablast)': 'megablast',
        'More dissimilar sequences (discontiguous megablast)': 'dc-megablast',
        'Somewhat similar sequences (blastn)': 'blastn'
    }
    task = task_mapping.get(task, task)

    print('{}: Starting {} for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), task, filename))

    db_folder = Path(blastn_database).joinpath('db')

    # Create a folder for subset BLASTN results
    blastn_subset_folder = project_folder.joinpath('subsets')
    os.makedirs(blastn_subset_folder, exist_ok=True)

    # Run BLASTN in parallel across all subsets
    print('{}: Your database: {}'.format(datetime.datetime.now().strftime('%H:%M:%S'), Path(blastn_database).stem))
    Parallel(n_jobs=n_cores, backend='threading')(delayed(blastn_parallel)(fasta_file, n_subsets, blastn_subset_folder, blastn_exe, db_folder, i, print_lock, task, max_target_seqs, masking) for i, fasta_file in enumerate(fasta_files))

    # Write log file with database and task information
    with open(project_folder.joinpath('log.txt'), 'w') as f:
            f.write(f'Blast executable:   {blastn_exe}\n')
            f.write(f'Blast database:     {blastn_database}\n')
            f.write(f'Output folder:      {project_folder}\n')
            f.write(f'Task:               {task}\n')
            f.write(f'Subset size:        {subset_size}\n')
            f.write(f'Max. target seq:    {max_target_seqs}\n')
            f.write(f'Masking:            {masking}\n')

    # Write OTU report
    with open(project_folder.joinpath('IDs.txt'), 'w') as f:
        for record in SeqIO.parse(query_fasta, "fasta"):
            f.write(record.id + '\n')

    print('{}: Finished {} for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), task, filename))

    # Merge BLASTN results (CSV files) into a single Parquet file with Snappy compression
    csv_files = glob.glob('{}/*.csv'.format(str(blastn_subset_folder)))
    col_names = ['unique ID', 'Sequence ID', 'Similarity', 'evalue']
    csv_files = [i for i in csv_files if os.path.getsize(i) != 0]
    df = pd.concat((pd.read_csv(f, header=None, sep=';;', names=col_names, engine='python').fillna('NAN') for f in csv_files))
    df['Sequence ID'] = df['Sequence ID'].astype(str)

    table = pa.Table.from_pandas(df)
    pq.write_table(table, accession_table, compression='snappy')

    # Remove temporary subset fasta folder
    shutil.rmtree(subset_folder)

if __name__ == '__main__':
    main()



