# apscale
Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data

# apscale-blast

## Introduction
Apscale-blast is part of the metabarcoding pipeline [apscale](https://github.com/DominikBuchner/apscale/tree/main) and is used
for the taxonomic assignment of OTUs or ESVs.

Programs used:
* [blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html) (blastn module)

Input:
* sequence fasta files (containing OTUs/ESVs) generated with Apscale (or any other metabarcoding pipeline)

Output:
* taxonomy table, log files

## Installation

Apscale-blast can be installed on all common operating systems (Windows, Linux, MacOS).
Apscale-blast requires Python 3.10 or higher and can be easily installed via pip in any command line:

`pip install apscale_blast`

To update apscale-blast run:

`pip install --upgrade apscale_blast`

## Distributables

Apscale-blast is also available as executable file that does not require the installation of python.
Please make sure to add the executable to the whitelist of your anti-virus software, otherwise the apscale-blast
will be blocked or deleted.

Distributables will soon be available [here](https://seafile.rlp.net/d/754a25aa76e44b2381b6/).

### Further dependencies - blast+

Apscale-blast calls blast+ (blastn) in the background. It should be installed and be in PATH to be executed
from anywhere on the system. However, a PATH to the executable can also be provided within the command.

You can find the latest blast+ executables and further information on the installation [here](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).

IMPORTANT: Please cite the blast+ software.

Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: Architecture and applications. BMC Bioinformatics, 10, 421. https://doi.org/10.1186/1471-2105-10-421

### Available databases

Apscale-blast uses pre-compiled databases. These databases will be tested and should prevent user-error. However, custom databases can also be created using the <i> database_creator.py </i> script (soon available here).

The pre-compiled databases are available under the following [server](https://seafile.rlp.net/d/754a25aa76e44b2381b6/) and will be updated regularly.

IMPORTANT: Please cite the used database accordingly!

#### Midori2

Leray, M., Knowlton, N., & Machida, R. J. (2022). MIDORI2: A collection of quality controlled, preformatted, and regularly updated reference databases for taxonomic assignment of eukaryotic mitochondrial sequences. Environmental DNA, 4(4), 894–907. https://doi.org/10.1002/edn3.303

#### Unite

Nilsson, R. H., Larsson, K.-H., Taylor, A. F. S., Bengtsson-Palme, J., Jeppesen, T. S., Schigel, D., Kennedy, P., Picard, K., Glöckner, F. O., Tedersoo, L., Saar, I., Kõljalg, U., & Abarenkov, K. (2019). The UNITE database for molecular identification of fungi: Handling dark taxa and parallel taxonomic classifications. Nucleic Acids Research, 47(D1), Article D1. https://doi.org/10.1093/nar/gky1022

When using the fungi database, please cite it as follows:

Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2024): UNITE general FASTA release for eukaryotes 2. Version 04.04.2024. UNITE Community. https://doi.org/10.15156/BIO/2959335

Includes global and 3% distance singletons.

When using the all eukaryotes database, please cite it as follows:

Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2024): UNITE general FASTA release for eukaryotes 2. Version 04.04.2024. UNITE Community. https://doi.org/10.15156/BIO/2959335

Includes global and 3% distance singletons.

When using the fungi database, please cite it as follows:

Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2024): UNITE general FASTA release for Fungi 2. Version 04.04.2024. UNITE Community. https://doi.org/10.15156/BIO/2959333

Includes global and 3% distance singletons.

#### SILVA

Quast, C., Pruesse, E., Yilmaz, P., Gerken, J., Schweer, T., Yarza, P., Peplies, J., & Glöckner, F. O. (2013). The SILVA ribosomal RNA gene database project: Improved data processing and web-based tools. Nucleic Acids Research, 41(D1), D590–D596. https://doi.org/10.1093/nar/gks1219

#### pr2

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L., Boutte, C., Burgaud, G., de Vargas, C., Decelle, J., del Campo, J., Dolan, J. R., Dunthorn, M., Edvardsen, B., Holzmann, M., Kooistra, W. H. C. F., Lara, E., Le Bescot, N., Logares, R., … Christen, R. (2013). The Protist Ribosomal Reference database (PR2): A catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with curated taxonomy. Nucleic Acids Research, 41(Database issue), D597–D604. https://doi.org/10.1093/nar/gks1160

#### diat.barcode

Rimet, F., Gusev, E., Kahlert, M., Kelly, M. G., Kulikovskiy, M., Maltsev, Y., Mann, D. G., Pfannkuchen, M., Trobajo, R., Vasselon, V., Zimmermann, J., & Bouchez, A. (2019). Diat.barcode, an open-access curated barcode library for diatoms. Scientific Reports, 9(1), Article 1. https://doi.org/10.1038/s41598-019-51500-6



## How to use

#### Pip installation

* Open a terminal, powershell, or command line window. Then, execute the following command:

`apscale_blast blastn`

#### Executable

* Open a terminal, powershell, or command line window. Navigate to the folder where the executable is located. Then, execute the following command:

Windows: `./apscale_blast.exe blastn`

MacOS: `./apscale_blast blastn`

#### Blastn search

* You will be asked to provide the FULL PATH to your database. Enter the PATH in the command line. For example:

`/Users/tillmacher/Downloads/MIDORI2_UNIQ_NUC_GB260_srRNA_BLAST`

* Next, you will be asked to provide the FULL PATH to your query fasta file that contains the OTU/ESV sequences. Enter the PATH in the command line. For example:

`/Users/tillmacher/Downloads/test_data_OTUs.fasta`

* Apscale-blast will split your fasta into smaller subsets [DEFAULT: 100].
  
* Each subset fasta file will be seperately blasted in parallel [DEFAULT: CPU count - 1].
  
* A new folder is created containing the unfiltered blastn results:

<pre>
├───test_data_OTUs.fasta
└───test_data_OTUs
    ├───IDs.txt
    ├───log.txt
    ├───test_data_OTUs.parquet.snappy
    └───subsets
        ├───subset_1_blastn.csv
        ├───subset_2_blastn.csv
        ├───subset_3_blastn.csv
        └───subset_4_blastn.csv
</pre>

#### Blastn results filtering

* Execute apscale-blast again and call the filter module:

`apscale_blast filter`

or 

Windows: `./apscale_blast.exe filter`

MacOS: `./apscale_blast filter`

* You will be asked to provide the FULL PATH to your database. Enter the PATH in the command line. For example:

`/Users/tillmacher/Downloads/MIDORI2_UNIQ_NUC_GB260_srRNA_BLAST`

* Next, you will be asked to provide the FULL PATH to the output folder that contains the unfiltered blastn results. Enter the PATH in the command line. For example:

`/Users/tillmacher/Downloads/test_data_OTUs`

* Apscale-blast will now filter the blastn results according to the following criteria:
  
- By e-value (the e-value is the number of expected hits of similar quality which could be found just by chance):
- The hit(s) with the lowest e-value are kept (the lower the e-value the better).
- By taxonomy:
- Hits with the same taxonomy are dereplicated.
- Hits are adjusted according to thresholds (default: species >=98%, genus >=95%, family >=90%, order >=85%) and dereplicated.
- Hits with still conflicting taxonomy are set back to the most recent common taxonomy
- OTU without matches in the blastn search are re-added as 'No Match'

* The results are collected and written to an new Excel file.

<pre>
├───test_data_OTUs.fasta
└───test_data_OTUs
    ├───IDs.txt
    ├───log.txt
    ├───test_data_OTUs_taxonomy.xlsx  <-- Filtered taxonomic identification results!
    ├───test_data_OTUs.parquet.snappy
    └───subsets
        ├───subset_1_blastn_filtered.xlsx
        ├───subset_1_blastn.csv
        ├───subset_2_blastn_filtered.xlsx
        ├───subset_2_blastn.csv
        ├───subset_3_blastn_filtered.xlsx
        ├───subset_3_blastn.csv
        ├───subset_4_blastn_filtered.xlsx
        └───subset_4_blastn.csv
</pre>

## Options

### apscale_blast blastn

options:

  -h, --help -> show this help message and exit
  
  -database DATABASE -> PATH to blastn database.
  
  -apscale_gui APSCALE_GUI -> Can be ignored: Only required for APSCALE-GUI.
  
  -blastn_exe BLASTN_EXE -> PATH to blast executable. [DEFAULT: blastn]
  
  -query_fasta QUERY_FASTA -> PATH to fasta file.
  
  -n_cores N_CORES -> Number of cores to use. [DEFAULT: CPU count - 1]
  
  -task TASK -> Blastn task: blastn, megablast, or dc-megablast. [DEFAULT: blastn]
  
  -out OUT -> PATH to output directory. A new folder will be created here. [DEFAULT: ./]
  
  -subset_size SUBSET_SIZE -> Number of sequences for each subset of the query fasta. [DEFAULT: 100]
                        
  -max_target_seqs MAX_TARGET_SEQS -> Number of hits retained from the blast search. Larger numbers will increase runtimes and required storage space [DEFAULT: 20]


### apscale_blast filter

options:

  -h, --help -> show this help message and exit
  
  -database DATABASE -> PATH to blastn database.
  
  -apscale_gui APSCALE_GUI -> Can be ignored: Only required for APSCALE-GUI.
  
  -blastn_folder BLASTN_FOLDER -> PATH to blastn folder for filtering.
  
  -thresholds THRESHOLDS -> Taxonomy filter thresholds. [DEFAULT: 97,95,90,87,85]
  
  -n_cores N_CORES -> Number of cores to use. [DEFAULT: CPU count - 1]










