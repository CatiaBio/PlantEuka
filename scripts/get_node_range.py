import subprocess
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument('name', help= "name list")
args = parser.parse_args()

f = open(args.name, "r")
f = f.read().splitlines()
name_list = list(f)


clade_info = {}
for name in name_list:
    out = subprocess.check_output(f'/home/ctools/vg/bin/vg stats -r /home/projects/animalia_mito/data/{name}/circ_{name}.vg', shell=True)
    ran = out.decode('utf-8')
    r_n = ran.split('\t')[1:]

    for i in r_n: 
        nodes = i.split(':')
    
    node_range = [int(x) for x in nodes]

    
    clade_info[name] = node_range



print(clade_info)




    





