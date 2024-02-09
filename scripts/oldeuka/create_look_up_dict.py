import argparse
import pickle
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fams', help = """list of available families created with ls -I "*.txt" > ava_family.txt""")
args = parser.parse_args()

f = open(args.fams, "r")
f = f.read().splitlines()
fams = list(f)


id_dict = defaultdict(list)
for i in fams:
    ass = []

    with open(f"/home/projects/animalia_mito/data/{i}/{i}_lin.txt") as file:
        for j in file:
            x = j.split(":")[0]
            if x not in ass:
                ass.append(x)
            else: 
                pass

    names = []        
    for l in range(len(ass)):
        names.append(f"{i}")


    for k,v in zip(ass,names): 
        id_dict[k].append(v)




with open('/home/projects/animalia_mito/scripts/look_up.pk', 'wb') as dict_file:
    pickle.dump(id_dict, dict_file)

