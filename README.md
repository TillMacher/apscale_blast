# APSCALE-blast

[![Downloads](https://static.pepy.tech/badge/apscale-blast)](https://pepy.tech/project/apscale-blast)

Taxonomic assignment of OTUs/ESVs with BLAST+, part of the metabarcoding pipeline [APSCALE](https://github.com/DominikBuchner/apscale).

* **Input:** a FASTA file with OTU/ESV sequences (from APSCALE or any other pipeline) and a reference database
* **Output:** a taxonomy table (`.xlsx`) with the best hit per sequence, a flag for ambiguous assignments and the thresholds used

## Installation

The easiest way is the conda environment of the [APSCALE installer](https://github.com/TillMacher/apscale_installer), which installs APSCALE-blast and all its dependencies.

Alternatively, install via pip (Python ≥ 3.10):

```
pip install apscale_blast
```

With pip, [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) must be installed separately and be available in your PATH.

## Usage

```
# Taxonomic assignment
apscale_blast -db ./MIDORI2_srRNA_BLAST -q ./ESVs.fasta -o ./my_run

# Re-BLAST sequences without a good match against a second database
apscale_blast -db ./MIDORI2_srRNA_BLAST -db2 ./my_local_db -q ./ESVs.fasta -o ./my_run

# Re-BLAST an existing result (results are saved next to the table as *_reblast.xlsx)
apscale_blast -reblast_only ./my_run/my_run_taxonomy.xlsx -db2 ./my_local_db -q ./ESVs.fasta

# Validate species assignments with GBIF occurrences (e.g. Germany, Luxembourg, France)
apscale_blast gbif -countries DE,LU,FR -q ./my_run/my_run_taxonomy.xlsx

# Only replace synonyms by their accepted names (GBIF backbone)
apscale_blast gbif -synonyms -q ./my_run/my_run_taxonomy.xlsx

# Calculate database-specific similarity thresholds
apscale_blast dboptimiser -db ./MIDORI2_srRNA_BLAST

# Create your own database (graphical interface)
apscale_blast dbcreator
```

All options are listed with `apscale_blast -h` (or e.g. `apscale_blast gbif -h`).

### How it works

1. The sequences are BLASTed (blastn) against the database in parallel subsets.
2. The best hits are filtered by similarity, e-value and reference rating.
3. Each assignment is trimmed to the taxonomic level supported by its similarity.
4. If several taxa remain, they are resolved with flags: F1 dominant species, F2 two species of one genus, F3 several species of one genus, F4 trimming to the most recent common ancestor.

### Similarity thresholds

* `apscale_blast dboptimiser` calculates thresholds per taxonomic level that best separate the similarity to the same taxon from the similarity to the closest other taxon (leave-one-out on the database). Use `-fp_weight 2` (or higher) for more conservative thresholds.
* If the database contains a `db_thresholds.xlsx` (created with `apscale_blast dboptimiser`), its database-wide thresholds are used.
* With `-threshold_level class` (or order, family, ...), each sequence gets the thresholds of its class. Classes based on fewer than 50 sequences use the database-wide thresholds.
* With `-thresholds 98,95,90,87,85` (species, genus, family, order, class), your own values are used.
* Without `db_thresholds.xlsx`, the default `98,95,90,87,85` is used.

## Databases

Pre-compiled databases are [available here](https://seafile.rlp.net/d/c172d076de1e4c45b594/) and are updated regularly. Custom databases can be created with `apscale_blast dbcreator`. Each database folder contains a `citation.txt`.

## Citation

Please cite BLAST+ and the database(s) you used:

* **BLAST+:** Camacho et al. (2009) BMC Bioinformatics 10: 421. https://doi.org/10.1186/1471-2105-10-421
* **MIDORI2:** Leray et al. (2022) Environmental DNA 4: 894–907. https://doi.org/10.1002/edn3.303
* **UNITE:** Abarenkov et al. (2024) Nucleic Acids Research 52: D791–D797. https://doi.org/10.1093/nar/gkad1039 (and the DOI of the [UNITE release](https://unite.ut.ee/repository.php) you used)
* **SILVA:** Quast et al. (2013) Nucleic Acids Research 41: D590–D596. https://doi.org/10.1093/nar/gks1219
* **PR2:** Guillou et al. (2013) Nucleic Acids Research 41: D597–D604. https://doi.org/10.1093/nar/gks1160
* **Diat.barcode:** Rimet et al. (2019) Scientific Reports 9: 15116. https://doi.org/10.1038/s41598-019-51500-6
* **ITS2 vascular plants:** Quaresma et al. (2024) Scientific Data 11: 129. https://doi.org/10.1038/s41597-024-02962-5
* **trnL (CRUX / Anacapa):** Curd et al. (2019) Methods in Ecology and Evolution 10: 1469–1475. https://doi.org/10.1111/2041-210X.13214
* **BOLD Systems:** Ratnasingham et al. (2024) BOLD v4, in DNA Barcoding: Methods and Protocols, pp. 403–441, plus the DOI of the [BOLD data package](https://boldsystems.org/about/citation/) you used
