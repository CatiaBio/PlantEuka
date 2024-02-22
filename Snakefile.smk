# Define paths to Python scripts
scripts = {
    "taxonomy": "scripts/get_taxonomy.py",
    "lineage": "scripts/get_lineage.py",
    "download": "scripts/download_fasta_v2.py",
    "organize": "scripts/organize_fasta.py",
    "taxID": "scripts/get_id_taxID.py",
    "update_mapping": "scripts/update_mapping_list.py"
}

rule all:
    input: 
        "data/chloroplast/organized",
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
        """
        python {scripts[taxonomy]} \
            --nodes_file {input.nodes} \
            --names_file {input.names} \
            --output_file {output} \
        """

rule generate_lineage:
    input: 
        "data/taxonomy.tsv"
    output:
        "data/lineage.tsv"
    shell: 
        """
        python {scripts[lineage]} \
            --taxonomy_file {input} \ 
            --output_file {output} \
        """

rule download_chloroplast_fasta:
    output:
        mapping="data/mapping_id_species_cp.txt",
        directory=directory("data/chloroplast/raw")  
    params:
        query="plants[filter] AND refseq[filter] AND chloroplast[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download]} \
            --search_str '{params.query}' \
            --mapping_file {output.mapping} \
            --output_dir {output.directory} \
        """ 

rule download_mitochondrion_fasta:
    output:
        mapping="data/mapping_id_species_mt.txt",
        directory=directory("data/mitochondrion/raw")  
    params:
        query="plants[filter] AND refseq[filter] AND mitochondrion[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download]} \
            --search_str '{params.query}' \
            --mapping_file {output.mapping} \
            --output_dir {output.directory} \
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
            --input_folder {input.folder_path}
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
            --input_folder {input.folder_path}
            --output_file {output.updated_mapping}
        """

rule check_contaminations:
    input:
        updated_mapping_mt = "data/combined_sequences_mt.mapping",
        updated_mapping_cp = "data/combined_sequences_cp.mapping",
        combined_sequences_mt = "data/combined_sequences_mt.fa.gz",
        combined_sequences_cp = "data/combined_sequences_cp.fa.gz"
    shell: 
        "conterminator dna {input.combined_sequences_mt} {input.updated_mapping_mt} combined_sequences_mt tmp_mt"
        "conterminator dna {input.combined_sequences_cp} {input.updated_mapping_cp} combined_sequences_cp tmp_cp"


