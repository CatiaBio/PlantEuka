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


# Define paths for tools
tools = {"stretcher": "/home/ctools/EMBOSS-6.6.0/emboss/stretcher",
    "mafft": "/home/ctools/Mafft/bin/mafft",
    "fasttree": "/home/ctools/FastTree/FastTree",
    "vg": "/home/ctools/vg_1.44.0/bin/vg"
}

rule all:
    input: 
        "data/combined_sequences_cp.mapping",
        "data/combined_sequences_mt.mapping"

rule download_taxdump_accession2taxid:
    output: 
        directory=directory("other")
    shell: 
        """ 
        wget -bqc -P other https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz 
        wget -bqc -P other https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
        """
    
rule download_cp_genomes:
    params:
        organelle = "chloroplast"
    shell: 
        """
        python {scripts[download_genomes]} \
        {params.organelle} 
        """ 

rule download_mt_genomes:
    params:
        organelle = "mitochondrion"
    shell: 
        """
        python {scripts[download_genomes]} \
        {params.organelle} 
        """ 

rule generate_lineage_taxonomy:
    shell:
        """
        python {scripts[taxonomy_lineage]} 
        """ 

rule generate_accession_taxid:
    shell:
        """
        python {scripts[accession_taxid]} 
        """ 

     
rule sort_mt_genomes:
    params:
        organelle = "mitochondrion"
    shell: 
        """
        python {scripts[sort]} \
        {params.organelle} 
        """

rule sort_cp_genomes:
    params:
        organelle = "chloroplast"
    shell:
    """
        python {scripts[sort]} \
        {params.organelle} 
        """

rule clean_genomes:
    input:
    output:
    shell: 
        """
        python {scripts[clean]} 
        """

