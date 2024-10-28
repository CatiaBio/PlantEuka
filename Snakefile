rule all:
    input: 
        "other/taxonomy.tsv",
        "other/lineage.tsv",
        "other/accession_taxid.txt"

rule download_taxdump_accession2taxid:
    output: 
        "other/taxdump.tar.gz",
        "other/nucl_gb.accession2taxid.gz"
    shell: 
        """ 
        wget -bqc -P other https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 
        wget -bqc -P other https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        """

rule unpack_taxdump:
    input:
        "other/taxdump.tar.gz"
    output:
        "other/nodes.dmp",
        "other/names.dmp"
    shell:
        """
        tar -zxvf {input} -C other nodes.dmp names.dmp
        """

rule generate_lineage_taxonomy:
    input:
        "other/nodes.dmp",
        "other/names.dmp"
    output:
        "other/taxonomy.tsv",
        "other/lineage.tsv"
    script:
        "scripts/taxonomy_lineage.py"

rule get_accession_list:
    input:
        "ncbi_info.txt"
    output:
        "other/accessions.txt"
    script:
        "scripts/get_accession_list.py"

rule download_genomes:
    input:
        ncbi_info="ncbi_info.txt",
        accession_list="other/accessions.txt"
    output:
        genomes_dir=directory("genomes/original"),  # Directory where genomes will be saved
        updated_downloaded_list="downloaded_accessions_list.txt"  # Updated file after downloading new accessions
    params:
        genomes_dir="genomes/original"  # Specify the directory for genomes explicitly if needed
    script:
        "scripts/download_genomes_.py"

rule generate_accession_taxid:
    input:
        "other/nucl_gb.accession2taxid.gz",
        "other/accessions.txt"
    output:
        "other/accession_taxid.txt"
    shell:
        "scripts/accession_taxid.sh"

     
# rule sort_cp_genomes:
#     params:
#         organelle = "chloroplast"
#     shell:
#     """
#         python {scripts[sort]} \
#         {params.organelle} 
#         """

# rule clean_genomes:
#     input:
#     output:
#     shell: 
#         """
#         python {scripts[clean]} 
#         """

