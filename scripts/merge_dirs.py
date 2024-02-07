import os
import argparse
import sys
import gzip
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fams', help = """Path to list of available families created with ls -I "*.txt" > ava_family.txt""")
parser.add_argument('-n', '--number', help= "add the number of sequences that should be the minimum for a family")
parser.add_argument('-c', '--create', help= "Takes an int, if the int is equal to 1 new family directories will be created, any other number and that step is passed")
args = parser.parse_args()

f = open(args.fams, "r")
f = f.read().splitlines()
fams = list(f)

wrong = []

for i in fams:
    with open(f"/home/projects/animalia_mito/data/{i}/{i}_lin.txt", "r") as lin:
                data = lin.readline()

                if f"{i}" in data: 
                    data = data.split("'")
                    
                    try:
                        fam = data.index(f"{i}")
                        
                        fam_name = data[fam - 2]
                        #print(fam_name)
                        #print(len(fam_name))
                        if len(fam_name == 1):
                            wrong.append(f"{i}")
                        
                        #if fam_name == "clade" or "group" or "RTA":
                        #    fam_name = data[fam - 4]
                        #else:
                         #   pass
                    except:
                        wrong.append(f"{i}")
#print(fam_name)

def check_families(fams, number_of_seqs):
    too_small = []
    os.system("touch /home/projects/animalia_mito/scripts/small_fams.txt")
    for i in fams:
        try:
    	    with gzip.open(f"/home/projects/animalia_mito/data/{i}/{i}_reps.fa.gz", "rt") as fam_file:
                seq = SeqIO.parse(fam_file, "fasta")
                n = 0
                for m in seq:
                    n += 1

                if n < number_of_seqs:
                    print(f"Genus {i} has not enough sequences.")
                    too_small.append(i)
        except:
            #pass
            with gzip.open(f"/home/projects/animalia_mito/data/{i}/{i}.fa.gz", "rt") as fam_file:
                seq = SeqIO.parse(fam_file, "fasta")
                n = 0
                for m in seq:
                    n += 1

                if n < number_of_seqs:
                    print(f"Family {i} has not enough sequences.")
                    too_small.append(i)
    return too_small


def sort_leftovers(update, new):
    os.system("touch /home/projects/animalia_mito/scripts/ava_next_fam.txt")
    for i in update:
        try:
            
            with open(f"/home/projects/animalia_mito/data/{i}/{i}_lin.txt", "r") as lin:
                data = lin.readline()

                if f"{i}" in data: 
                    data = data.split("'")
                    fam = data.index(f"{i}")
                    
                    fam_name = data[fam - 2]
                    #print(fam_name)
                        #if fam_name == "clade" or "group" or "RTA" or "lineage" or "sedis":
                            #print(fam_name)#fam_name = data[fam - 4]
                elif len(fam_name == 1):
                    pass

                if os.path.isdir(f"/home/projects/animalia_mito/data/{fam_name}/"):
                    os.system(f"""mv /home/projects/animalia_mito/data/{i} /home/projects/animalia_mito/data/{fam_name}""")
                    print(f"Directory already exists. {i} is moved to {fam_name}", file = sys.stdout)
                else:
                    print("Directory not found. Waiting for next round.", file = sys.stdout)

                if new == 1:
                    if os.path.isdir(f"/home/projects/animalia_mito/data/{fam_name}/"):
                        os.system(f"""mv /home/projects/animalia_mito/data/{i} /home/projects/animalia_mito/data/{fam_name}""")
                        print(f"Directory already exists. {i} is moved to {fam_name}", file = sys.stdout)
                    else:
                        os.system(f"""mkdir /home/projects/animalia_mito/data/{fam_name}""")
                        os.system(f"""mv /home/projects/animalia_mito/data/{i} /home/projects/animalia_mito/data/{fam_name}/""")
                        os.system(f"""echo {fam_name} >> ava_family.txt""")
                        print(f"New directory {fam_name} has been created. {i} has been moved into the new family directory", file = sys.stout)
                else:
                    pass
        except:
            print(f"Could not find lineage file for {i}.")


                

def create_fam_fasta():
    os.system("ls /home/projects/animalia_mito/data/ > /home/projects/animalia_mito/scripts/all.txt")
    f = open("/home/projects/animalia_mito/scripts/all.txt", "r")
    f = f.read().splitlines()
    f = list(f)
    for i in f:
        if os.path.isfile(f"/home/projects/animalia_mito/data/{i}/{i}.fa.gz"):
            pass
        else:
            os.system(f"zcat /home/projects/animalia_mito/data/{i}/*/*_reps.fa.gz /home/projects/animalia_mito/data/{i}/*/*/*_reps.fa.gz /home/projects/animalia_mito/data/{i}/*/*/*/*_reps.fa.gz /home/projects/animalia_mito/data/{i}/*/*/*/*/*_reps.fa.gz > /home/projects/animalia_mito/data/{i}/{i}.fa")
            os.system(f"gzip /home/projects/animalia_mito/data/{i}/{i}.fa")
            print(f"New multi-fasta file {i}.fa.gz was created.", file = sys.stdout)

def create_fam_lin():
    os.system("ls /home/projects/animalia_mito/data/ > /home/projects/animalia_mito/scripts/all.txt")
    f = open("/home/projects/animalia_mito/scripts/all.txt", "r")
    f = f.read().splitlines()
    f = list(f)
    for i in f:
        if os.path.isfile(f"/home/projects/animalia_mito/data/{i}/{i}_lin.txt"):
            pass
        else:
            os.system(f"cat /home/projects/animalia_mito/data/{i}/*/*_lin.txt /home/projects/animalia_mito/data/{i}/*/*/*_lin.txt /home/projects/animalia_mito/data/{i}/*/*/*/*_lin.txt /home/projects/animalia_mito/data/{i}/*/*/*/*/*_lin.txt  > /home/projects/animalia_mito/data/{i}/{i}_lin.txt")
            print(f"New lineage file {i}_lin.txt was created.", file = sys.stdout)



print(wrong)
#print(len(wrong))
#update = check_families(fams, int(args.number))
#round1 = sort_leftovers(update, int(args.create))
new = create_fam_fasta()
lin = create_fam_lin()
