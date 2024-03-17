# Define paths to Python scripts
scripts = {
    "taxonomy": "scripts/generate_taxonomy_list.py",
    "lineage": "scripts/generate_lineage_list.py",
    "download": "scripts/download_genomes.py",
    "organize": "scripts/organize_fasta.py",
    "taxID": "scripts/generate_id_taxid.py",
    "clean": "scripts/clean_genomes.sh"
}

rule all:
    input: 
        "data/combined_sequences_cp.mapping",
        "data/combined_sequences_mt.mapping"

rule download_taxdump:
    output: temp("taxdump.tar.gz")
    shell: 
        "wget -O {output} https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

rule extract_taxdump:
    input: "taxdump.tar.gz"
    output: 
        nodes=temp("data/taxdump/nodes.dmp"),
        names=temp("data/taxdump/names.dmp")
    shell:
        """
        mkdir -p other
        tar -zxvf {input} -C other nodes.dmp names.dmp
        """

rule generate_taxonomy:
    input:
        nodes="other/nodes.dmp",
        names="othernames.dmp"
    output: 
        "other/taxonomy.tsv"
    shell: 
        """
        python {scripts[taxonomy]} \
            {input.nodes} \
            {input.names} \
            {output} \
        """

rule generate_lineage:
    input: 
        "other/taxonomy.tsv"
    output:
        "other/lineage.tsv"
    shell: 
        """
        python {scripts[lineage]} \
            {input} \ 
            {output} \
        """

rule download_chloroplast_genomes:
    output:
        id_species_list="genomes/chloroplast/id_species_cp.tsv",
        directory=directory("genomes/chloroplast/unsorted")  
    params:
        query="plants[filter] AND refseq[filter] AND chloroplast[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download]} \
            {params.query} \
            {output.mapping} \
            {output.directory} \
        """ 

rule download_mitochondrion_genomes:
    output:
        id_species_list="genomes/mitochondrion/id_species_mt.tsv",
        directory=directory("genomes/mitochondrion/unsorted")  
    params:
        query="plants[filter] AND refseq[filter] AND mitochondrion[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download]} \
            {params.query} \
            {output.mapping} \
            {output.directory} \
        """ 
     
rule organize_mitochondrion_genomes:
    input: 
        mapping="genomes/mitochondrion/id_species_mt.tsv",
        raw_directory="genomes/mitochondrion/unsorted",  
        lineage="other/lineage.tsv"
    output:
        directory=directory("genomes/mitochondrion/sorted")  
    params:
        input_dir="genomes/mitochondrion/unsorted",  
        output_dir="genomes/mitochondrion/sorted"  
    shell: 
        """
        python {scripts[organize]} \
            {input.lineage} \
            {input.mapping} \
            {params.input_dir} \
            {params.output_dir}
        """

rule organize_chloroplast_genomes:
    input: 
        mapping="genomes/chloroplast/id_species_mt.tsv",
        raw_directory="genomes/chloroplast/unsorted",  
        lineage="data/lineage.tsv"
    output:
        directory=directory("genomes/chloroplast/sorted")  
    params:
        input_dir="genomes/chloroplast/unsorted",  
        output_dir="genomes/chloroplast/sorted" 
    shell: 
        """
        python {scripts[organize]} \
            {input.lineage} \
            {input.mapping} \
            {params.input_dir} \
            {params.output_dir}
        """

rule get_chloroplast_taxID:
    input:
        mapping="genomes/chloroplast/id_species_cp.tsv", 
        taxonomy="other/taxonomy.tsv"
    output:
        mappingtaxID="genomes/id_taxid_cp.txt"
    shell: 
        """
        python {scripts[taxID]} \
            {input.mapping} \
            {input.taxonomy} \
            {output.mappingtaxID}
        """

rule get_mitochondrion_taxID:
    input:
        mapping="genomes/mitochondrion/id_species_mt.tsv", 
        taxonomy="other/taxonomy.tsv"
    output:
        mappingtaxID="genomes/mitochondrion/id_taxID_mt.tsv"
    shell: 
        """
        python {scripts[taxID]} \
            {input.mapping} \
            {input.taxonomy} \
            {output.mappingtaxID}
        """

rule clean_chloroplast_genomes: 
    input:
        genomes= "genomes/chloroplast/sorted",
    output:
        directory=directory("logs")
    params:
        log_clean = "logs/genome_cleanup_cp.log"
    shell: 
        """
        ./{scripts[clean]} \
            {input.genomes} \
            {params.log_clean} \   
        """

rule clean_mitochondrion_genomes:  
    input:
        genomes= "genomes/mitochondrion/sorted/",
    output:
        directory=directory("logs")
    params:
        log_clean = "logs/genome_cleanup_mt.log"
    shell: 
        """
        ./{scripts[clean]} \
            {input.genomes} \
            {params.log_clean} \   
        """

# rule check_contaminations_mt:
#     input:
#         updated_mapping_mt = "data/combined_sequences_mt.mapping",
#         combined_sequences_mt = "data/combined_sequences_mt.fa.gz"
#     params:
#         output_dir=directory("data/tmp_mt")
#     shell: 
#         "tools/conterminator/conterminator dna {input.combined_sequences_mt} {input.updated_mapping_mt} combined_sequences_mt {output.output_dir}"

# rule check_contaminations_cp:
#     input:
#         updated_mapping_cp = "data/combined_sequences_cp.mapping",
#         combined_sequences_cp = "data/combined_sequences_cp.fa.gz"
#     params:
#         output_dir=directory("data/tmp_cp")
#     shell: 
#         "tools/conterminator/conterminator dna {input.combined_sequences_cp} {input.updated_mapping_cp} combined_sequences_cp {output.output_dir}"

