import argparse
import multiprocessing
import os
import sys
import datetime
from pathlib import Path
from apscale_blast.a_blastn import main as a_blastn
from apscale_blast.b_filter import main as b_filter
from ete3 import NCBITaxa, Tree

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

def main():
    """
    APSCALE BLASTn suite
    Command-line tool to run and filter BLASTn searches.
    """

    # Introductory message with usage examples
    message = """
    APSCALE blast command line tool - v1.2.2
    Example commands:
    $ apscale_blast -h
    $ apscale_blast -db ./MIDORI2_UNIQ_NUC_GB259_srRNA_BLAST -q ./12S_apscale_ESVs.fasta -f Vertebrata
    
    Remember to update your local ete3 NCBI taxonomy regularly, if using the "remote" blastn!
    This can be performed by running:
    $ apscale_blast -u
    """
    print(message)

    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='APSCALE blast v1.2.2')

    # === Main settings ===
    main_settings = parser.add_argument_group("Main settings")
    main_settings.add_argument('-database', '-db', type=str, required=False, help='PATH to local database. Use "remote" to blast against GenBank (slow).')
    main_settings.add_argument('-query_fasta', '-q', type=str, help='PATH to fasta file.')
    main_settings.add_argument('-out', '-o', type=str, default='./', help='PATH to output directory. A new folder will be created here. [DEFAULT: ./]')

    # === General settings ===
    general_settings = parser.add_argument_group("General settings")
    general_settings.add_argument('-n_cores', type=int, default=multiprocessing.cpu_count() - 1, help='Number of CPU cores to use. [DEFAULT: CPU count - 1]')
    general_settings.add_argument('-task', type=str, default='blastn', help='Blastn task: blastn, megablast, or dc-megablast. [DEFAULT: blastn]')
    general_settings.add_argument('-subset_size', type=int, default=100, help='Number of sequences per query fasta subset. [DEFAULT: 100]')
    general_settings.add_argument('-max_target_seqs', type=int, default=20, help='Number of hits retained from the blast search. [DEFAULT: 20]')
    general_settings.add_argument('-thresholds', type=str, default='97,95,90,87,85', help='Taxonomy filter thresholds. [DEFAULT: 97,95,90,87,85]')

    # === Remote blast settings ===
    remote_blast = parser.add_argument_group("Remote blast settings")
    remote_blast.add_argument('-update_taxids', '-u', action='store_true', help='Update NCBI taxid backbone.')
    remote_blast.add_argument('-organism_filter', '-f', type=str, help='Comma-separated list of taxids or full names for remote blast filtering (e.g., "Mammalia,Actinopteri").')
    remote_blast.add_argument('-include_uncultured', action='store_true', help='Include uncultured/environmental sample sequences in the remote blast [DEFAULT=True].')

    # === Advanced settings ===
    advanced_settings = parser.add_argument_group("Advanced settings")
    advanced_settings.add_argument('-blastn_exe', type=str, default='blastn', help='PATH to blast executable. [DEFAULT: blastn]')
    advanced_settings.add_argument('-masking', action='store_false', help='Activate masking. [DEFAULT=False]')
    advanced_settings.add_argument('-disable_headless', action='store_false', help='Show/hide the chromium browser [DEFAULT=True]')

    # Parse the arguments
    args = parser.parse_args()

    # Handle taxonomy update
    if args.update_taxids:
        print("Updating NCBI taxonomy database...")
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()
        print("Taxonomy database updated successfully.")
        if Path('./taxdump.tar.gz').exists():
            os.remove('./taxdump.tar.gz')
            print('Removed taxdmup.tar.gz')
        return  # Exit after updating taxonomy

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
    project_folder = args.out  # Use the output directory specified by the user

    # Handle the 'blastn' command
    if args.disable_headless == False:
        headless = False
    else:
        headless = True
    continue_blast = False

    # Convert db to Path
    database = Path(args.database.strip('"'))
    if str(database) == 'remote':
        database = 'remote'

    if ' ' in database:
        print('\nError: Your database PATH contains white spaces! Please move your database to a different folder!')

    if database == 'remote' and not args.organism_filter:
        print('\nError: Remote blast requires at least one organism to filter for!')
        print('Example: -f Mammalia, Actinopteri')
        return

    organism_mask = []
    if database == 'remote':
        for organism in args.organism_filter.replace(' ', '').split(','):
            organism_mask.append(organism_filter(organism))

        if len(organism_mask) != 0:
            print(f'{datetime.datetime.now().strftime("%H:%M:%S")}: Filtering remote blast for {len(organism_mask)} organism(s):')
            for mask in organism_mask:
                if 'not found!' in mask:
                    print(f'  Error: Please check your input: {mask}\n')
                    return
                else:
                    print(f'  {mask}')

    if args.query_fasta:
        # Run the BLASTn function
        continue_blast = a_blastn(args.blastn_exe,
                 args.query_fasta.strip('"'),
                 database,
                 project_folder,
                 args.n_cores,
                 args.task,
                 args.subset_size,
                 args.max_target_seqs,
                 args.masking,
                 headless,
                 organism_mask,
                 args.include_uncultured
                                  )
    else:
        print('\nError: Please provide a fasta file!')

    # Handle the 'filter' command
    if continue_blast == False:
        print('\nNot all fasta subsets have been processed yet!')
    elif not os.path.isfile(Path(project_folder).joinpath('log.txt')):
        print('\nError: Could not find the BLAST results folder!')
    else:
        # Run the filter function
        b_filter(project_folder, database, args.thresholds, args.n_cores)

# Run the main function if script is called directly
if __name__ == "__main__":
    main()