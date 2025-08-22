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

# Lock for thread-safe print statements
print_lock = threading.Lock()

# Ensure multiprocessing works when frozen (e.g. for executables)
multiprocessing.freeze_support()

def fasta_subset(fasta_file, subset_size):
    """
    Splits a large fasta file into smaller subsets for parallel processing.

    Args:
        fasta_file (str): Path to the input fasta file.
        subset_size (int): Number of sequences per subset file.

    Returns:
        Path: Path to the directory containing the subsets.
    """
    print('{}: Creating subset(s) from fasta file.'.format(datetime.datetime.now().strftime('%H:%M:%S')))

    subset_size = int(subset_size)
    fasta_file = Path(fasta_file)

    # Create a new directory for subsets
    subset_folder = Path(fasta_file.parent).joinpath('fasta_subsets')
    os.makedirs(subset_folder, exist_ok=True)

    # Delete existing subset files, if any
    for f in glob.glob(str(subset_folder / '*.fasta')):
        os.remove(f)

    chunk_fasta_files = []
    i, n = 1, 1

    # Splitting fasta file into chunks
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            chunk_fasta = '{}/subset_{}.fasta'.format(subset_folder, i)
            if chunk_fasta not in chunk_fasta_files:
                chunk_fasta_files.append(chunk_fasta)

            with open(chunk_fasta, 'a') as output_handle:
                SeqIO.write(record, output_handle, 'fasta')

            # Create new chunk after reaching subset_size
            if n == subset_size:
                n = 1
                i += 1
            else:
                n += 1

    print('{}: Created {} subset(s) from fasta file.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i))
    return subset_folder

def accession2taxonomy(df_1, taxid_dict, col_names_2, db_name):
    """
    Maps accession numbers to taxonomy using a dictionary.

    Args:
        df_1 (pd.DataFrame): DataFrame containing accession data.
        taxid_dict (dict): Dictionary mapping accessions to taxonomy.
        col_names_2 (list): Column names for the output DataFrame.
        db_name (str): Name of the database.

    Returns:
        pd.DataFrame: DataFrame with taxonomy information appended.
    """
    df_2_list = []
    for row in df_1.values.tolist():
        ID_name, accession = row[0], row[1]
        evalue, similarity = row[-1], row[-2]

        taxonomy = taxid_dict.get(accession, ['NoMatch'] * 7)
        df_2_list.append([ID_name] + taxonomy + [similarity, evalue])

    df_2 = pd.DataFrame(df_2_list, columns=col_names_2)
    return df_2

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

def remote_blast(fasta_file, n_subsets, blastn_subset_folder, blastn_exe, db_folder, i, print_lock, task, max_target_seqs, masking, headless, tmp_folder, organism_mask, include_uncultured):

    print('{}: Starting remote blast for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))

    # Define output files
    blastn_json = blastn_subset_folder.joinpath(Path(fasta_file).stem + '_' + task + '.json')
    blastn_csv = blastn_subset_folder.joinpath(Path(fasta_file).stem + '_' + task + '.csv')

    # Create tmp file
    t = datetime.datetime.now().strftime('%D_%H_%M_%S').replace('/', '_')
    tmp_file = tmp_folder.joinpath(f'{t}_search.txt')

    # Combine all sequences from the FASTA file into a single query string
    query = ""
    for record in SeqIO.parse(fasta_file, "fasta"):
        query += f">{record.id}\n{str(record.seq).upper()}\n"

    # Skip if output already exists
    if os.path.isfile(blastn_csv) and os.path.getsize(blastn_csv) > 0:
        with print_lock:
            print('{}: Skipping {} (already exists and is not empty).'.format(datetime.datetime.now().strftime('%H:%M:%S'),
                                                             blastn_csv.stem))
    else:

        # headless = 'False'
        # query = 'ACGT'
        # blastn_subset_folder = '/Users/tillmacher/Desktop/APSCALE_projects/test_apscale/8_esv_table/quatsch_tax'
        # i = 1

        with sync_playwright() as p:
            # Launch the browser
            browser = p.chromium.launch(headless=headless)
            context = browser.new_context()

            # Create a new page
            page = context.new_page()

            # Set up a directory to save the downloaded file
            download_dir = blastn_subset_folder
            os.makedirs(download_dir, exist_ok=True)

            # Define query details for naming
            query_name = f'Query_{i}'  # Extracted from input or hardcoded
            timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")  # Current date and time

            # Navigate to the BLAST page
            url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch#"
            page.goto(url)

            # 1) Paste the query sequence into the input textarea
            textarea_selector = "#seq"
            page.fill(textarea_selector, query)
            time.sleep(random.randrange(1,5))

            # 1.1) Paste the organism to include in the input area
            c = 0
            for mask in organism_mask:
                if c == 0:
                    textarea_selector = '#qorganism'
                    page.fill(textarea_selector, mask)
                    c += 1
                else:
                    page.wait_for_selector("#addOrg", state="visible")  # Ensure it's visible
                    page.click("#addOrg")
                    textarea_selector = f'#qorganism{c}'
                    page.fill(textarea_selector, mask)
                    c+=1


            # 1.2) Select to exclude unculured samples sequences
            if include_uncultured == False:
                exclude_uncultured_selector = "label[for='exclSeqUncult']"
                page.click(exclude_uncultured_selector)

            time.sleep(random.randrange(1,5))

            # Step 2: Choose "discontiguous megablast"
            available_algorithms = {'megablast': 'megaBlast', 'dc-megablast': 'discoMegablast', 'blastn': 'blastn'}
            selected_algorithm = available_algorithms[task]
            # Interact with the label associated with the radio button
            label_selector = f"label[for='{selected_algorithm}']"  # Select the label for the radio button
            # Ensure the label is visible and click it
            label = page.locator(label_selector)
            label.scroll_into_view_if_needed()
            label.click()
            # Verify that the radio button is selected
            radio_button_selector = f"input#{selected_algorithm}"
            radio_button = page.locator(radio_button_selector)
            time.sleep(random.randrange(1,5))

            # 3) Run BLAST
            blast_button_selector = "#blastButton1 > input.blastbutton"
            page.click(blast_button_selector)

            # 4) Wait for the page to load the results
            page.wait_for_selector("#allDownload", state="attached", timeout=900 * 1000)  # Wait for the download section to appear
            time.sleep(random.randrange(5,12))  # Wait for BLAST results to be processed

            # 5) Click the "Download All" button to reveal options
            download_all_button_selector = "#ulDnldAl"
            page.click(download_all_button_selector)

            # 6) Wait for the menu to open
            page.wait_for_selector("#allDownload[aria-hidden='false']", state="visible")  # Ensure the menu is now visible

            # 7) Click the "Single-file JSON" download link
            single_file_json_selector = "a.xgl[href*='FORMAT_TYPE=JSON2_S']"

            # Expect download and click the link
            with page.expect_download() as download_info:
                page.click(single_file_json_selector)

            # Get the downloaded file path
            download = download_info.value
            download_path = download.path()  # Path to the downloaded file

            # Construct a meaningful filename
            blastn_json_path = os.path.join(download_dir, blastn_json)

            if os.path.isfile(download_path):
                # Rename and move the file
                if not download_path.exists():
                    print('Error: Could not find download!')
                    print(download_path)
                    print(blastn_json_path)
                shutil.move(download_path, blastn_json_path)

                # Allow some time for the file to move before closing the browser
                time.sleep(random.randrange(5,15))

                with open(blastn_json_path, "r") as file:
                    json_data = json.load(file)

                # Convert JSON to csv table
                json_to_csv(json_data, blastn_json_path)

                # Close the browser
                browser.close()

                # Sleep long to reduce over-stressing the server
                time.sleep(random.randrange(25, 35))

                print('{}: Finished remote blast for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))

                tmp_file.touch()

            else:
                # Close the browser
                browser.close()

                # Sleep long to reduce over-stressing the server
                time.sleep(random.randrange(5, 15))

                print('{}: Failed remote blast for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1, n_subsets))

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
    elif masking == "No":
        # Run the BLASTN command
        subprocess.call([blastn_exe, '-task', task, '-db', str(db_folder), '-query', str(fasta_file),
                         '-num_threads', str(1), '-max_target_seqs', str(max_target_seqs),
                         '-dust', 'no', '-soft_masking', 'false',
                         '-outfmt', '6 delim=;; qseqid sseqid pident evalue', '-out', str(blastn_csv)])
        with print_lock:
            print('{}: Finished blastn for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1,
                                                                 n_subsets))

    else:
        # Run the BLASTN command
        subprocess.call([blastn_exe, '-task', task, '-db', str(db_folder), '-query', str(fasta_file),
                         '-num_threads', str(1), '-max_target_seqs', str(max_target_seqs),
                         '-outfmt', '6 delim=;; qseqid sseqid pident evalue', '-out', str(blastn_csv)])
        with print_lock:
            print('{}: Finished blastn for subset {}/{}.'.format(datetime.datetime.now().strftime('%H:%M:%S'), i + 1,
                                                                 n_subsets))

def main(blastn_exe, query_fasta, blastn_database, project_folder, n_cores, task, subset_size, max_target_seqs, masking, headless, organism_mask, include_uncultured):
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

    filename = Path(query_fasta).stem.replace('.', '_').replace(' ', '_')

    print('{}: Starting {} for \'{}\''.format(datetime.datetime.now().strftime('%H:%M:%S'), task, filename))
    db_folder = Path(blastn_database).joinpath('db')

    # Create a folder for subset BLASTN results
    blastn_subset_folder = project_folder.joinpath('subsets')
    os.makedirs(blastn_subset_folder, exist_ok=True)

    continue_blast = True
    do_not_ask_again = False

    if blastn_database == 'remote':
        # Run remote blast NOT IN PARALLEL!
        limit = 10
        potential_sequences = limit * subset_size
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Using remote BLAST.')
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Please note that requests will be rate-limited to {limit} requests per day ({potential_sequences} sequences) to avoid overloading the server.')
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Subsets will be processed sequentially to minimize server strain.')
        print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: For larger datasets (1000+ query sequences), please use a local database.')

        # Define output directory for tmp files
        script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
        tmp_folder = script_dir.joinpath('tmp')
        os.makedirs(tmp_folder, exist_ok=True)

        for i, fasta_file in enumerate(fasta_files):
            # Check how many runs were performed today
            today = datetime.datetime.now().strftime('%D').replace('/', '_')
            n_runs = len(glob.glob(str(tmp_folder.joinpath(f'{today}*.txt'))))

            if n_runs > limit and not do_not_ask_again:
                print('')
                print(f"{datetime.datetime.now().strftime('%H:%M:%S')}: More than {n_runs} remote blasts have been requested today!")
                print("  It is advised to continue blasting tomorrow to prevent penalties on the IP address!")
                print("  Tip: Use one of the local databases for more than 2000 ESVs.")
                print("  The NCBI online blast module is free for everybody and should be used fairly.")

                print("  Options:")
                print("  1 = Continue (will ask again for each file)")
                print("  2 = Do not ask again and continue for all files (NOT RECOMMENDED)")
                print("  3 = Stop, relax, and continue tomorrow (RECOMMENDED)")

                answer = input("Enter your choice: ").strip()

                if answer == '1':  # Continue only for this run
                    continue_blast = True
                    do_not_ask_again = False
                elif answer == '2':  # Do not ask again
                    continue_blast = True
                    do_not_ask_again = True
                else:  # Stop the process
                    return False

            if continue_blast == True or do_not_ask_again == True:
                remote_blast(
                    fasta_file,
                    n_subsets,
                    blastn_subset_folder,
                    blastn_exe,
                    db_folder,
                    i,
                    print_lock,
                    task,
                    max_target_seqs,
                    masking,
                    headless,
                    tmp_folder,
                    organism_mask,
                    include_uncultured
                )
                print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Conducted {n_runs+1}/{limit} remote blasts today.')

    else:
        print('{}: Your database: {}'.format(datetime.datetime.now().strftime('%H:%M:%S'), Path(blastn_database).stem))
        # Run BLASTN in parallel across all subsets
        Parallel(n_jobs=n_cores, backend='threading')(delayed(blastn_parallel)(
            fasta_file, n_subsets, blastn_subset_folder, blastn_exe, db_folder, i, print_lock, task, max_target_seqs, masking
        ) for i, fasta_file in enumerate(fasta_files))

    # Write log file with database and task information
    with open(project_folder.joinpath('log.txt'), 'w') as f:
            f.write(f'Blast executable:   {blastn_exe}\n')
            f.write(f'Blast database:     {blastn_database}\n')
            f.write(f'Output folder:      {project_folder}\n')
            f.write(f'Task:               {task}\n')
            f.write(f'Subset size:        {subset_size}\n')
            f.write(f'Max. target seq:    {max_target_seqs}\n')
            f.write(f'Masking:            {masking}\n')
            f.write(f'Organism mask:      {organism_mask}\n')
            f.write(f'Exclude uncultured: {include_uncultured}\n')

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
    pq.write_table(table, project_folder.joinpath('{}.parquet.snappy'.format(filename)), compression='snappy')

    # Remove temporary subset fasta folder
    shutil.rmtree(subset_folder)

    # tell b_filter to continue or not
    return continue_blast

if __name__ == '__main__':
    main()
