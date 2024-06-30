# PlantEuka 🌱🌿🌻

**PlantEuka** is an extension of **Euka**, developed by Nicola et al. (https://doi.org/10.1111/2041-210X.14214), which uses a taxon-based pangenome graph for analyzing ancient environmental DNA (aeDNA). While **Euka** has shown great potential, it is limited to tetrapod and arthropod mitochondrial genomes (mitogenomes).

The **PlantEuka** project aims to develop a customized plant database module for **Euka**. This enhancement will improve the accuracy and representation of plant species, enabling the identification of plant species using a taxon-based pangenome graph.

The primary goal of this project is to enhance our understanding of plant biodiversity and their interactions in various environments. These improvements will allow **Euka** to provide more reliable data for ecological and biological research.

## How to Use the PlantEka

Download the repository and follow the steps. Ensure all scripts are executed from the main directory of PlantEuka.

### download_genomes.py

This script interfaces with NCBI’s GenBank database to search and retrieve sequences matching a specified query string.

**Usage:** 
```bash
scripts/download_genomes.py <organelle>
```

### Taxonomy files 

Download taxonomy files, taxdump.tar.gz and nucl_gb.accession2taxid.gz from NCBI's FTP repository 

```bash
wget -bqc -P other https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
wget -bqc -P other https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
```

### accession_taxid.sh

This script matches downloaded accession numbers with their corresponding taxon IDs (taxids) using files from the NCBI database.

**Usage:** 
```bash
scripts/accession_taxid.sh
```

### taxonomy_lineage.py

This script extracts taxonomy and lineage information from `nodes.dmp` and `names.dmp` files, creating `taxonomy.tsv` and `lineage.tsv`.

**Usage:** 
```bash
scripts/taxonomy_lineage.py
```

### sort_genomes.py

This script organizes genomic data into directories by taxonomic rank (genus, family, order), creating lists of sorted and unsorted accession numbers.

**Usage:** 
```bash
scripts/sort_genomes.py <organelle>
```

### clean_genomes.py

The script scans genomic sequences for non-standard characters and replaces them with a standard ’N’ for unknown nucleotides

**Usage:** 
```bash
scripts/clean_genomes.py <organelle>
```

Note: Still under development. More details will be added soon.
