import sys
import os
import argparse
from random import randint
from time import sleep
from download_green_algae_test import download_w_taxa

parser = argparse.ArgumentParser()
parser.add_argument('--path', help="Path to file list")
parser.add_argument('--search', help="Which search term: Chlorophyta?")
args = parser.parse_args()

genus_list = open(args.path, "r").read().splitlines()

for genus in genus_list:
    if os.path.isdir(f"/home/projects/MAAG/PlantEuka/cp_DB/data/{genus}/"):
        pass
    else:
        # Define search term for NCBI based on the argument provided
        search = ""
        if args.search == "Chlorophyta":
            search = f"""(("{genus}"[Organism] OR {genus}[All Fields]) AND chloroplast[All Fields] AND complete genome[All Fields]) AND refseq[filter]"""

        if search:  # Proceed only if the search term is defined
            dwt = download_w_taxa(search, genus)

            # Check if the search term returns a result; if not, the script skips to the next genus
            os.system(f"touch /home/projects/MAAG/PlantEuka/cp_DB/scripts/downloaded_{args.search}_genus.txt")
            try:
                check_1 = dwt.create_fasta()
                print("is this working?")
                print(check_1)

                sleep(randint(1,20))
                if check_1 == 0:
                    print(f"Genus {genus} did not return any records within NCBI's database for the specified search requirements.", file=sys.stdout)
                else:
                    dwt.sort_fasta()
                    dwt.exclude_samples()
                    try:
                        test = dwt.make_lineage_file()
                        if test == 0:
                            os.system(f"rm -r /home/projects/MAAG/PlantEuka/cp_DB/data/{genus}")
                            print(f"Wrongly added genus {genus} was removed", file=sys.stdout)
                        else:
                            os.system(f"gzip /home/projects/MAAG/PlantEuka/cp_DB/data/{genus}/*.fa")
                            os.system(f"echo {genus} >> /home/projects/MAAG/PlantEuka/cp_DB/scripts/downloaded_{args.search}_genus.txt")
                            print(f"Directory {genus} was created.", file=sys.stdout)
                    except:
                        print(f"Lineage file couldn't be created - check {genus}", file=sys.stdout)
            except:
                print(f"Genus {genus} did not return any records within NCBI's database for the specified search requirements.", file=sys.stdout)
