# Define paths to Python scripts
scripts = {
    "taxonomy": "scripts/generate_taxonomy_list.py",
    "lineage": "scripts/generate_lineage_list.py",
    "download_genomes": "scripts/download_genomes.py",
    "organize": "scripts/organize_genomes.py",
    "taxID": "scripts/generate_id_taxid.py",
    "clean": "scripts/clean_genomes.py",
    "merge": "scripts/merge_genomes.sh",
    "gene_info": "scripts/download_gene.py",
    "gene_rank": "scripts/match_gene_with_rank.py",
    "gene_start": "scripts/move_gene_to_start.py"
    "pairs_pairwise": "scripts/create_list_pairs_full_paths.sh"
    "pairs_pairwise_results": "scripts/write_pairs_parallel.sh"
}


# Define paths for tools
tools = {
    "stretcher" = "/home/ctools/EMBOSS-6.6.0/emboss/stretcher",
    "mafft" = "/home/ctools/Mafft/bin/mafft",
    "fasttree" = "/home/ctools/FastTree/FastTree"
      
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
        nodes=temp("other/taxdump/nodes.dmp"),
        names=temp("other/taxdump/names.dmp")
    shell:
        """
        mkdir -p other 
        tar -zxvf {input} -C other nodes.dmp names.dmp
        """

rule generate_taxonomy:
    input:
        nodes="other/nodes.dmp",
        names="other/names.dmp",
        plant_taxid = "33090"
    output: 
        "other/taxonomy.tsv"
    shell: 
        """
        python {scripts[taxonomy]} \
            {input.nodes} \
            {input.names} \
            {input.plant_taxid} \
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
        directory=directory("genomes/chloroplast/original")  
    params:
        query="plants[filter] AND refseq[filter] AND chloroplast[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download_genomes]} \
            {params.query} \
            {output.mapping} \
            {output.directory} \
        """ 

rule download_mitochondrion_genomes:
    output:
        id_species_list="genomes/mitochondrion/id_species_mt.tsv",
        directory=directory("genomes/mitochondrion/original")  
    params:
        query="plants[filter] AND refseq[filter] AND mitochondrion[filter] AND complete genome[Title]"
    shell: 
        """
        python {scripts[download_genomes]} \
            {params.query} \
            {output.mapping} \
            {output.directory} \
        """ 
     
rule organize_mitochondrion_genomes:
    input: 
        mapping="genomes/mitochondrion/id_species_mt.tsv",
        raw_directory="genomes/mitochondrion/original",  
        lineage="other/lineage.tsv"
    output:
        directory=directory("genomes/mitochondrion/sorted")  
    params:
        input_dir="genomes/mitochondrion/original",  
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
        raw_directory="genomes/chloroplast/original",  
        lineage="data/lineage.tsv"
    output:
        directory=directory("genomes/chloroplast/sorted")  
    params:
        input_dir="genomes/chloroplast/original",  
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
        genomes= "genomes/chloroplast/sorted"
    output:
        directory=directory("logs")
    params:
        log_clean = "logs/genome_cleanup_cp.log"
    shell: 
        """
        python {scripts[clean]} \
            {input.genomes} \
            {params.log_clean}   
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
        python \
        {scripts[clean]} \
        {input.genomes} \
        {params.log_clean}  
        """
################################################################################################################################################################################################
####### This part is an example of how to run the cluster in parallel (possibly will be removed)
### This was needed to run the program for the pairwise alignment in parallel
## Firstly a list of all the possible pairs is created (rule create_all_pairs_pairwise)
## Secondly a lits of all results names is created and added to the list. This will be used to name the pairwise result file containing the name of the pairs that were used. 

rule create_all_pairs_pairwise: 
# Creates a list with the paths to all possible pairs of fasta files within each group in a rank (genus,family,order)
    input:
        genomes_cp = "genomes/chloroplast/sorted"
        genomes_mt = "genomes/mitochondrion/sorted"
    output:
        directory=directory("other/pw_lists/chloroplast")
        directory=directory("other/pw_lists/mitochondrion")
    shell: 
        """
        shell \
        {scripts[pairs_pairwise]} \
        {input.genomes_cp} \
        {input.genomes_mt}   
        """

rule pairs_pairwise_results: 
# Adds to the list a path to the pairwise result to be used as output after running the pairwise program
    input:
        pairwise_cp = "other/pw_lists/chloroplast"
        pairwise_mt = "other/pw_lists/mitochondrion"
    shell: 
        """
        shell \
        {scripts[pairs_pairwise]} \
        {input.genomes_cp} \
        {input.genomes_mt}   
        """

rule run_pairwise_cp_parallel: 
# This is an example to run the list of pairs using stretcher in parallel
# list_serv_2_14 is a txt file with a list of the servers to use for parallel 
    input:
        pairwise_cp = "other/pw_lists/combined_parallel_cp.txt"
        pairwise_mt = "other/pw_lists/combined_parallel_mt.txt"
    shell: 
        """
        cat combined_parallel_cp.txt combined_parallel_mt.txt > combined_parallel.txt | \
        cat combined_parallel.txt | \
        awk '{print "nice -19 {tools[strectcher]} \
        -asequence <(zcat "$1") \
        -bsequence <(zcat "$2") \
        -gapopen 16 -gapextend 4 \
        -outfile >(gzip > "$3".gz)"}' | \
        parallel --slf list_serv_2_14  
        """
################################################################################################################################################################################################

# Example to run stretcher for pairwise 
rule run_pairwise_cp: 
    input:
        pair1 = "genomes/chloroplast/sorted/genus/Abies/NC_026892.1.fasta.gz"
        pair2 =  "genomes/chloroplast/sorted/genus/Abies/NC_057315.1.fasta.gz"
    output:
        pairwise_result = "results/pairwise/genus_Abies_NC_026892.1_NC_057315.1.stretcher.gz"
    shell:
        """"
        {tools[strectcher]} \
        -asequence <(zcat {input.pair1}) \
        -bsequence <(zcat {input.pair1}) \
        -gapopen 16 -gapextend 4 \
        -outfile >(gzip > {output.pairwise_result})
        """

# Example to run mafft 
    input:
        merged_genome = "genomes/chloroplast/sorted/merged/Abies/NC_026892.1.fasta.gz"
    output:
        msa = "results/msa/genus_Abies_NC_026892.1.fasta.gz"
    shell:
        """"
        {tools[mafft]} \
        --auto --thread -1 <(zcat {input.merged_genome}) \
        > >(gzip > {output.msa} )
        """

# Example to run fasttree 
rule run_mafft_cp: 
    input:
        msa = "results/msa/genus_Abies_NC_026892.1.fasta.gz"
    output:
        tree = "results/msa/genus_Abies_NC_026892.1.tree"
    shell:
        """"
        {tools[fasttree]} -nt -gtr -gamma \
        {input.msa} \
        > {output.tree}
        """