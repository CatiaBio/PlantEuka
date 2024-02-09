import os 
import argparse
import sys 


p = argparse.ArgumentParser()
p.add_argument('genus_list')
args = p.parse_args()

g = open(args.genus_list, "r")
g = g.read().splitlines()
genus = list(g)

def genus_to_family (genus):
	os.system("touch /home/projects/animalia_mito/scripts/ava_family.txt")
   	for i in genus:
		with open(f"/home/projects/animalia_mito/data/{i}/{i}_lin.txt", "r") as lin:
			data = lin.readline()
			if f"{i}" in data: 
				data = data.split("'")
				fam = data.index(f"{i}")
				fam_name = data[fam - 2]
				if fam_name == "clade" or "group" or "RTA":
					fam_name = data[fam - 4]
				else:
					pass

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


fams = genus_to_family(genus)

