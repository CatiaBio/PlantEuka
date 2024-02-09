#!/usr/bin/python

import subprocess
import re
import os
import csv
import json
import pandas as pd
import numpy as np
from scipy import stats
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('name', help= "name list") #Give it the all_clades list that is used to create the combo graph
args = parser.parse_args()

f = open(args.name, "r")
f = f.read().splitlines()
name_list = list(f)

clade_chunks = csv.writer(open("/home/projects/animalia_mito/scripts/clade_chunks5.csv", "w"))
for name in name_list:
    print(name)
    if os.path.isfile(f"""/home/projects/animalia_mito/data/{name}/circ_{name}.vg"""):
        os.system(f"""/home/ctools/vg_1.44.0/bin/vg chunk -x /home/projects/animalia_mito/data/{name}/circ_{name}.vg -s 1500 -c 0""")

        os.system(f"""/home/ctools/vg_1.44.0/bin/vg paths -x /home/projects/animalia_mito/data/{name}/circ_{name}.vg -L > paths.txt""")
        p = open("paths.txt", "r")
        p = p.read().splitlines()
        p_list = list(p)
        clade_bin_list = []
        comp_list2 = []
        for j in p_list:
            os.system(f"""ls -r *{j}* > chunks.txt""")
            c = open("chunks.txt", "r")
            c = c.read().splitlines()
            c_list = list(c)
            bins = []
            comp_list = []
            for m in c_list:
                
                os.system(f"""/home/ctools/vg_1.44.0/bin/vg view {m} -j > {m}.json""")

                # = open(f"{m}.json")
                
                
                with open(f"{m}.json", 'r', encoding='utf-8') as t:
                    #import pdb; pdb.set_trace()
                    data = json.loads(t.readlines()[0])
                    seq = ''
                    for d in data['node']:
                        seq+=d['sequence']

                    

                    seq_s = []
                    n=1
                    for index in range(0,len(seq),n):
                        seq_s.append(seq[index:index+n])

                    seq_series = pd.Series(seq_s)
                    count = seq_series.value_counts()
                    seq_complexity = stats.entropy(count)
                    comp_list.append(seq_complexity)
                    



                    out = subprocess.check_output(f"""/home/ctools/vg_1.44.0/bin/vg stats -r {m}""", shell=True)
                    ran = out.decode('utf-8')
                    r_n = re.findall('(\d+)', ran)
                    bins.append(r_n)


            
            if len(bins) > 1:
                bin = [[int(j) for j in i] for i in bins]
                bin.sort()
                check = bin[0]
            
                del bin[0]
                if check in bin:
                    index = bin.index(check)
                    del bin[index]
                else:
                    pass

                # adding a new first bin as the original one just takes the entire range of the graph. 
                check_new = []
                check_new.append(check[0])
                check_new.append(bin[0][0])
                bin.insert(0, check_new)
                

                flat_b = list(np.concatenate(bin). flat)
                b = [int(i) for i in flat_b]

                
                clade_bin_list.append(b)

            
            comp_list2.append(comp_list)   


        
        bin_df = pd.DataFrame(clade_bin_list)

        com_df = pd.DataFrame(comp_list2)

                
        
        bin_range_max = bin_df.max(axis = 0, skipna = True)
        bin_range_min = bin_df.min(axis = 0, skipna = True)
        bin_complexity = com_df.median(axis = 0, skipna = True)
        min_values = list(bin_range_min)[::2]
        print(len(min_values))
        max_values = list(bin_range_max)[1::2]
        print(len(max_values))
        median_complexity = list(bin_complexity)
        print(len(median_complexity))
        fin_range = []
        for g,h,k in zip(min_values, max_values, median_complexity):

            fin_range.append(g)
            fin_range.append(h)
            fin_range.append(k)
        

        
        fin_range.insert(0, name)
        print(fin_range)
        clade_chunks.writerow(fin_range)

        os.system("rm chunk_*")
            
              

    else: 
        pass 
    
    
