import sys
import os
import numpy as np
import pandas as pd
import argparse
from random import randint
from time import sleep
from downloader import download_w_taxa


parser = argparse.ArgumentParser()
parser.add_argument('--path', help= "Path to file list")
parser.add_argument('--search', help= "Which search term: Arthro, Tetra, Human?")
args = parser.parse_args()

artho = open(args.path, "r")
arth_new = artho.read().splitlines()
artho_new = list(arth_new)
check_genera = []

for i in artho_new:

    genus = i
    if os.path.isdir(f"/home/projects/animalia_mito/data/{genus}/"):
        pass
    else:
        #search term for NCBI
        if args.search == "Arthro":
            search = f"""(("{genus}"[Organism] OR {genus}[All Fields]) AND "{genus}"[title] AND mitochondrion[All Fields] AND complete[All Fields] AND genome[All Fields]) AND refseq[filter]""" 
        elif args.search == "Tetra":
            search = f"""(("{genus}"[Organism] OR {genus}[All Fields]) AND "{genus}"[title] AND mitochondrion[All Fields] AND complete[All Fields] AND genome[All Fields]) AND refseq[filter] NOT isolate[title]""" 

#f"""(10000[SLEN] : 20000[SLEN])
 #       AND "{genus}"[Organism]
  #      AND mitochondrion[All Fields]
   #     AND complete genome[title]
    #    AND "{genus}"[title]
     #   NOT voucher[title]
      #  NOT strain[title]
       # NOT clone[title]
       # NOT haplotype[title]
       # NOT breed[title]"""
        elif args.search == "Human":
            search = """HM771189[All Fields] OR MF621128[All Fields] OR DQ304984[All Fields] OR MF437265[All Fields] OR MF696041[All Fields] OR DQ341063[All Fields] OR MK059556[All Fields] OR MT554224[All Fields] OR JQ702641[All Fields] OR KJ919975[All Fields] OR JQ705183[All Fields] OR JQ705480[All Fields] OR KY595668[All Fields] OR AY339507[All Fields] OR MT795654[All Fields] OR KT780370[All Fields] AND (10000[SLEN] : 18000[SLEN])"""
        

        dwt = download_w_taxa(search, genus)
        
        #check if the search term returns a result, if not the class skips to the next genus
        os.system(f"touch /home/projects/animalia_mito/scripts/downloaded_{args.search}_genus.txt")
        try:
            check_1 = dwt.create_fasta()
            print("is this working?")
            print(check_1)
        
            
            sleep(randint(1,20))
            if check_1 == 0:
                print(f"Genus {genus} did not return any records within NCBI's database for the specified search requirements.", file = sys.stdout)
            else:
                sort_fasta = dwt.sort_fasta()
                po = dwt.exclude_samples()
                try:
                    test = dwt.make_lineage_file()
                    if test == 0:
                        os.system(f"rm -r /home/projects/animalia_mito/data/{genus}")
                        print(f"Wrongly added genus {genus} was removed", file = sys.stdout)
                    else:
                        os.system(f"gzip /home/projects/animalia_mito/data/{genus}/*.fa")
                        os.system(f"echo {genus} >> /home/projects/animalia_mito/scripts/downloaded_arthro_genus.txt")
                        print(f"Directory {genus} was created.", file = sys.stdout)
                except:
                    check_genera.append(f"{genus}")
                    print(f"Lineage file couldn't be created - check {genus}", file = sys.stdout)
        except:
            print(f"Genus {genus} did not return any records within NCBI's database for the specified search requirements.", file = sys.stdout)
         

print(check_genera, file = sys.stdout)
