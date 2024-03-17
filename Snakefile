# Define paths to Python scripts
scripts = {
    "taxonomy": "scripts/generate_taxonomy_list.py",
    "lineage": "scripts/generate_lineage_list.py",
    "download": "scripts/download_genomes.py",
    "organize": "scripts/organize_fasta.py",
    "taxID": "scripts/get_id_taxID.py",
    "update_mapping": "scripts/update_mapping_list.py"
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
        id_species_list="genomes/chloroplast/unsigned/id_species_cp.tsv",
        directory=directory("genomes/chloroplast/unsigned")  
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
        id_species_list="genomes/mitochondrion/unsigned/id_species_mt.tsv",
        directory=directory("genomes/mitochondrion/unsigned")  
    params:
        query="plants[filter] AND refseq[filter] AND mitochondrion[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download]} \
            {params.query} \
            {output.mapping} \
            {output.directory} \
        """ 
     
rule organize_mitochondrion_fasta:
    input: 
        mapping="data/mapping_id_species_mt.txt",
        raw_directory="data/mitochondrion/raw",  
        lineage="data/lineage.tsv"
    output:
        directory=directory("data/mitochondrion/organized")  
    params:
        input_dir="data/mitochondrion/raw",  
        output_dir="data/mitochondrion/organized"  
    shell: 
        """
        python {scripts[organize]} \
            --lineage_file {input.lineage} \
            --mapping_file {input.mapping} \
            --input_dir {params.input_dir} \
            --output_dir {params.output_dir}
        """

rule organize_chloroplast_fasta:
    input: 
        mapping="data/mapping_id_species_cp.txt",
        raw_directory="data/chloroplast/raw",  
        lineage="data/lineage.tsv"
    output:
        directory=directory("data/chloroplast/organized")  
    params:
        input_dir="data/chloroplast/raw",  
        output_dir="data/chloroplast/organized" 
    shell: 
        """
        python {scripts[organize]} \
            --lineage_file {input.lineage} \
            --mapping_file {input.mapping} \
            --input_dir {params.input_dir} \
            --output_dir {params.output_dir}
        """

rule get_chloroplast_taxID:
    input:
        mapping="data/mapping_id_species_cp.txt", 
        taxonomy="data/taxonomy.tsv"
    output:
        mappingtaxID="data/mapping_id_taxID_cp.txt"
    shell: 
        """
        python {scripts[taxID]} \
            --mapping_file {input.mapping} \
            --taxonomy_file {input.taxonomy} \
            --output_file {output.mappingtaxID}
        """

rule get_mitochondrion_taxID:
    input:
        mapping="data/mapping_id_species_mt.txt", 
        taxonomy="data/taxonomy.tsv"
    output:
        mappingtaxID="data/mapping_id_taxID_mt.txt"
    shell: 
        """
        python {scripts[taxID]} \
            --mapping_file {input.mapping} \
            --taxonomy_file {input.taxonomy} \
            --output_file {output.mappingtaxID}
        """

rule update_mapping_list_mt: 
    input:
        mapping= "data/mapping_id_taxID_mt.txt",
        folder_path = "data/mitochondrion/raw"
    output:
        updated_mapping = "data/combined_sequences_mt.mapping"
    shell: 
        """
        python {scripts[update_mapping]} \
            --mapping_file {input.mapping} \
            --input_folder {input.folder_path} \
            --output_file {output.updated_mapping}
        """

rule update_mapping_list_cp: 
    input:
        mapping = "data/mapping_id_taxID_cp.txt",
        folder_path = "data/chloroplast/raw"
    output:
        updated_mapping = "data/combined_sequences_cp.mapping"
    shell: 
        """
        python {scripts[update_mapping]} \
            --mapping_file {input.mapping} \
            --input_folder {input.folder_path} \
            --output_file {output.updated_mapping} 
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

