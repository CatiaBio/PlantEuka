rule all:
    input: 
        "other/taxonomy.tsv",
        "other/lineage.tsv",
        "other/accession_taxid.txt",
        "chloroplast/other/sorted_accessions.txt",
        "chloroplast/genomes/sorted",
        "chloroplast/other/genome_cleanup.log",
        "chloroplast/other/all_pairs_parallel.txt",
        "chloroplast/results/pairwise/pairwise_results.tsv",
        "chloroplast/results/pairwise/filtered_results.tsv",
        "chloroplast/genomes/merged",
        "chloroplast/other/msa_paths.txt",
        "chloroplast/results/msa"

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

rule sort_cp_genomes:
    input:
        lineage="other/lineage.tsv",
        accession_taxid="other/accession_taxid.txt",
        genomes_dir=directory("genomes/original")
    output:
        sorted_list="chloroplast/other/sorted_accessions.txt",
        unsorted_list="chloroplast/other/unsorted_accessions.txt",
        sorted_dir=directory("chloroplast/genomes/sorted")
    params:
        organelle="chloroplast"
    script:
        "scripts/organize_genomes.py"

rule clean_genomes:
    input:
        sorted_list="chloroplast/other/sorted_accessions.txt",
        sorted_dir="chloroplast/genomes/sorted"
    output:
        log_file="chloroplast/other/genome_cleanup.log"
    params:
        organelle="chloroplast"
    script:
        "scripts/clean_genomes.py"

rule generate_pairwise_list:
    input:
        log_file="chloroplast/other/genome_cleanup.log",
        sorted_dir="chloroplast/genomes/sorted"
    output:
        pairs_list="chloroplast/other/all_pairs_parallel.txt"
    params:
        organelle="chloroplast"
    shell:
        "scripts/pairwise_list_parallel.sh {params.organelle}"

rule run_stretcher_alignments:
    input:
        pairs_list="chloroplast/other/all_pairs_parallel.txt"
    output:
        results_dir=directory("chloroplast/results/pairwise/stretcher")
    params:
        organelle="chloroplast"
    shell:
        """
        mkdir -p {output.results_dir}
        while IFS=' ' read -r seq1 seq2 output_file; do
            if [ -f "$seq1" ] && [ -f "$seq2" ]; then
                gunzip -c "$seq1" > temp_seq1.fasta
                gunzip -c "$seq2" > temp_seq2.fasta
                stretcher -asequence temp_seq1.fasta -bsequence temp_seq2.fasta -gapopen 16 -gapextend 4 -outfile "$output_file"
                rm temp_seq1.fasta temp_seq2.fasta
            fi
        done < {input.pairs_list}
        """

rule gather_pairwise_results:
    input:
        results_dir="chloroplast/results/pairwise/stretcher"
    output:
        results_file="chloroplast/results/pairwise/pairwise_results.tsv"
    params:
        organelle="chloroplast"
    script:
        "scripts/pairwise_results.py"

rule filter_pairwise_results:
    input:
        results_file="chloroplast/results/pairwise/pairwise_results.tsv"
    output:
        filtered_file="chloroplast/results/pairwise/filtered_results.tsv"
    params:
        organelle="chloroplast"
    script:
        "scripts/pairwise_filter.py"

rule merge_genomes:
    input:
        filtered_file="chloroplast/results/pairwise/filtered_results.tsv",
        sorted_dir="chloroplast/genomes/sorted"
    output:
        merged_dir=directory("chloroplast/genomes/merged")
    params:
        organelle="chloroplast"
    shell:
        """
        mkdir -p {output.merged_dir}
        scripts/merge_genomes.sh {input.sorted_dir} {output.merged_dir}
        """

rule generate_msa_paths:
    input:
        merged_dir="chloroplast/genomes/merged"
    output:
        msa_paths="chloroplast/other/msa_paths.txt"
    params:
        organelle="chloroplast"
    script:
        "scripts/msa_paths_parallel.py"

rule run_mafft_alignments:
    input:
        msa_paths="chloroplast/other/msa_paths.txt"
    output:
        results_dir=directory("chloroplast/results/msa")
    params:
        organelle="chloroplast"
    shell:
        """
        mkdir -p {output.results_dir}
        while read input_file output_file; do
            if [ -f "$input_file" ]; then
                mafft --auto --thread -1 "$input_file" > "$output_file"
            fi
        done < {input.msa_paths}
        """

