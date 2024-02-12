# Define paths to Python scripts
scripts = {
    "generate_taxonomy_list": "scripts/viridiplantae_taxonomy_list.py",
    "generate_genus_list": "scripts/viridiplantae_genus_list.py"
}

rule all:
    input: "data/viridiplantae_genus_list.tsv"

rule download_taxdump:
    output: "taxdump.tar.gz"
    shell: "wget -O {output} https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

rule extract_taxdump:
    input: "taxdump.tar.gz"
    output: 
        nodes="data/taxdump/nodes.dmp",
        names="data/taxdump/names.dmp"
    shell:
        """
        mkdir -p data/taxdump
        echo 'Extracting specific files from {input}...'
        tar -zxvf {input} -C data/taxdump nodes.dmp names.dmp
        echo 'Extraction of specific files complete.'
        rm taxdump.tar.gz
        touch {output.nodes} {output.names}
        """

rule generate_viridiplantae_taxonomy_list:
    input:
        nodes="data/taxdump/nodes.dmp",
        names="data/taxdump/names.dmp"
    output: "data/viridiplantae_taxonomy.tsv"
    shell: "python {scripts[generate_taxonomy_list]} --nodes_file {input.nodes} --names_file {input.names} --output_file {output}"

rule generate_viridiplantae_genus_list:
    input: 
        "data/viridiplantae_taxonomy.tsv"
    output:
        "data/viridiplantae_genus_list.tsv"
    shell: 
        "python {scripts[generate_genus_list]} --taxonomy_file {input} --output_file {output}"
