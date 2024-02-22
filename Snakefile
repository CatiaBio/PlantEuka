# Define paths to Python scripts
scripts = {
    "taxonomy": "scripts/get_taxonomy.py",
    "lineage": "scripts/get_lineage.py",
    "download": "scripts/download_fasta_v2.py",
    "organize": "scripts/organize_fasta.py",
    "taxID": "scripts/get_taxID.py"
}

# Update the 'all' rule to reflect the actual final outputs of the workflow
rule all:
    input: 
        #"data/chloroplast/organized",
        "data/mitochondrion/organized"

rule download_taxdump:
    output: temp("taxdump.tar.gz")
    shell: 
        "wget -O {output} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

rule extract_taxdump:
    input: "taxdump.tar.gz"
    output: 
        nodes=temp("data/taxdump/nodes.dmp"),
        names=temp("data/taxdump/names.dmp")
    shell:
        """
        mkdir -p data/taxdump
        tar -zxvf {input} -C data/taxdump nodes.dmp names.dmp
        """

rule generate_taxonomy:
    input:
        nodes="data/taxdump/nodes.dmp",
        names="data/taxdump/names.dmp"
    output: 
        "data/taxonomy.tsv"
    shell: 
        "python {scripts[taxonomy]} --nodes_file {input.nodes} --names_file {input.names} --output_file {output}"

rule generate_lineage:
    input: 
        "data/taxonomy.tsv"
    output:
        "data/lineage.tsv"
    shell: 
        "python {scripts[lineage]} --taxonomy_file {input} --output_file {output}"

# rule download_chloroplast_fasta:
#     input: 
#         query="plants[filter] AND refseq[filter] AND chloroplast[filter] AND complete genome[Title]"
#     output:
#         mapping="data/mapping_id_species_cp.txt",
#         directory="data/chloroplast"
#     shell: 
#         "python {scripts[download]} --search_str '{input.query}' --mapping_file {output.mapping} --output_dir {output.directory}"

rule download_mitochondrion_fasta:
    output:
        mapping="data/mapping_id_species_mt.txt",
        directory=directory("data/mitochondrion/raw")  # Adjusted path to avoid overlap
    params:
        query="plants[filter] AND refseq[filter] AND mitochondrion[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download]} \
            --search_str '{params.query}' \
            --mapping_file {output.mapping} \
            --output_dir {output.directory} \
        """ 

# rule organize_chloroplast_fasta:
#     input: 
#         mapping="data/mapping_id_species_cp.txt",
#         lineage="data/lineage.tsv"
#     output:
#         directory="data/chloroplast/organized"
#     shell: 
#         "python {scripts[organize]} --lineage_file {input.lineage} --mapping_file {input.mapping} --output_dir {output.directory}"


rule get_mitochondrion_taxID:
    input:
        mapping="data/mapping_id_species_mt.txt", 
        taxonomy="data/taxonomy.tsv"
    output:
        mappingtaxID="data/mapping_id_species_taxID_mt.txt"
    shell: 
        """
        python {scripts[taxID]} \
            --mapping_file {input.mapping} \
            --taxonomy_file {input.taxonomy} \
            --output_file {output.mappingtaxID}
        """

rule organize_mitochondrion_fasta:
    input: 
        mapping="data/mapping_id_species_mt.txt",
        raw_directory="data/mitochondrion/raw",  # Ensure this matches the directory structure used for downloading
        lineage="data/lineage.tsv"
    output:
        directory=directory("data/mitochondrion/organized")  # Using directory() to indicate the output is a directory
    params:
        input_dir="data/mitochondrion/raw",  # Specify the input directory for the FASTA files
        output_dir="data/mitochondrion/organized"  # Specify the output directory for organized FASTA files
    shell: 
        """
        python {scripts[organize]} \
            --lineage_file {input.lineage} \
            --mapping_file {input.mapping} \
            --input_dir {params.input_dir} \
            --output_dir {params.output_dir}
        """


