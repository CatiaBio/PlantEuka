#!/usr/bin/python

import pickle
import subprocess
import re
import numpy as np
import csv
import os
import argparse
import pdb


parser = argparse.ArgumentParser()
parser.add_argument('name', help= "name list")
args = parser.parse_args()

f = open(args.name, "r")
f = f.read().splitlines()
name_list = list(f)

# import the look up dictonary for the assession numbers

file = open("/home/projects/animalia_mito/scripts/look_up.pk", "rb")
id_dict = pickle.load(file)

print(id_dict)


head_tail = []
for name in name_list:
        out = subprocess.check_output(f'/home/ctools/vg/bin/vg stats -r /home/projects/animalia_mito/data/{name}/circ_{name}.vg', shell=True)
        ran = out.decode('utf-8')
        r_n = ran.split('\t')[1:]

        for i in r_n: 
                nodes = i.split(':')
        
        node_range = [int(x) for x in nodes]

        head_tail.append(node_range)

print(head_tail)




clade_info = csv.writer(open("/home/projects/animalia_mito/vgan_dev/share/euka_dir/clade_information.csv", "w"))



cc = 0  
for i in head_tail:
    print(i)
    clade = []
    no = subprocess.check_output(f"""/home/ctools/vg/bin/vg trace -x /home/projects/animalia_mito/vgan_dev/share/euka_dir/euka_db.vg -G /home/projects/animalia_mito/vgan_dev/share/euka_dir/euka_db.gbwt -n {i[0]} -a {i[1]} | /home/ctools/vg/bin/vg view /dev/stdin | grep "P"| head -1 | cut -f 2 """, shell=True)
    print(no)
    nono = no.decode('utf-8')
    #slice after numberic occurrence
    assession_no = nono.split('_')[0]
    if (len(assession_no) > 2 ):
        pass
    else:
        assession_no = "_".join(nono.split('_')[:2])
    print(assession_no)
    c_name = id_dict[f"{assession_no}"]
    print(c_name)
    print(type(c_name))
    for m in c_name:
        # if cc == 26:
        #     pdb.set_trace()
        print(len(c_name))
        print(m)
        dis = subprocess.check_output(f"/home/ctools/opt/R-4.0.5/bin/Rscript get_pairwise_dist.R {m}", shell=True)
        print("Hello I pass this")
        print(dis)
        av = dis.decode('utf-8')
        av_dis = av.split(' ')[1]
        av_dis = av_dis.replace("\n", "")
        
        clade.append(cc)
        clade.append(f"{m}")
        clade.append(av_dis)

        print(clade)
        clade_info.writerow(clade)
        

      
    cc += 1
        


#os.system("sed 's/,/ /g' /home/projects/animalia_mito/vg_res/clade_information.csv > /home/projects/animalia_mito/vg_res/clade_info.csv")


#print(clade_info)
#clade.append(cc)
        #clade.append(f"{m}")

#clade_dump = np.array(clade_info)
#print(clade_dump)
#np.savetxt(sys.stdout.buffer, clade_dump)
#clade_dump.tofile("/home/projects/animalia_mito/vg_res/clade_information.csv", sep="\t")

#with open("/home/projects/animalia_mito/vg_res/clade_info.csv", "w+") as file:
    #file.write(str(clade_dump))
    #file.close()



    

