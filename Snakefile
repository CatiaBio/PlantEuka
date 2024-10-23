# Define paths to Python scripts
scripts = {
    "accession_taxid": "scripts/accession_taxid.sh",
    "taxonomy_lineage": "scripts/taxonomy_lineage.py",
    "download_genomes": "scripts/download_genomes.py",
    "sort": "scripts/sort_genomes.py",
    "clean": "scripts/clean_genomes.py",
    "pw_list": "pw_list_parallel.sh",
    "merge": "scripts/merge_genomes.sh"
}


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
    shell:
        """
        python {scripts[taxonomy_lineage]} 
        """ 

rule download_cp_genomes:
    params:
        email = "catiacarmobatista@gmail.com",
        api_key = "09a80c0d55826098773b6a9e63c0514f5508"
    output:
        "other/accessions.txt"
    shell: 
        """
        python {scripts[download_genomes]} {params.email} {params.api_key}        
        """ 

rule generate_accession_taxid:
    input:
        "other/nucl_gb.accession2taxid.gz",
        "other/accessions.txt"
    output:
        "other/accession_taxid.txt"
    shell:
        """
        python {scripts[accession_taxid]} 
        """ 

     
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

