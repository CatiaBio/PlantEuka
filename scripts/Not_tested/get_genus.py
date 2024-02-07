import pandas as pd

taxa = pd.read_csv("../data/taxa.tab", sep="\t", header = None, error_bad_lines=False, names=['a', 'b', 'name', 'c', 'rank', 'taxid'])

path = "/home/biokate/thesis_unix/test_data/scripts/"

chlorophyta_genus = []
for index, row in taxa.iterrows():
    if row['rank'] == 'genus' and '3041' in row['taxid']:
        chlorophyta_genus.append(row['name'])
        
print(len(chlorophyta_genus))
chloro = open((path + "chloro_genus"), "w")
for e in chloro:
    chloro.write(e + "\n")
chloro.close()