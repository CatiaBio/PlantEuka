from Bio import Entrez,SeqIO
import os
from ete3 import NCBITaxa


class download_w_taxa:
    """Downloads all complete mitogenomes from NCBI for Athropods and Tetrapods.
    Download batches at genus level, all fasta files will be accompanied by a 
    taxon file. The storage order will be based on the phylogeny of the taxa."""
    Entrez.email = "n.alexandra.vogel@gmail.com"
    Entrez.api_key = "f815551132f17a6be4e57546520bd8b46007"

    # Holds the download location of NCBI's taxonomic database
    ncbi = NCBITaxa()
    #update function if necessary
    #ncbi.update_taxonomy_database()

    path = "/home/biokate/thesis_unix/test_data/scripts"

    def __init__(self, search_statement, genus):
        self.search_statement = search_statement
        self.genus = genus
        self.path = f"/home/projects/animalia_mito/data/{self.genus}/{self.genus}.fa"
        self.path2 = f"/home/projects/animalia_mito/data/{self.genus}/{self.genus}_reps.fa"
     

    # creates fasta file based off of the search statement given 
    # if fasta list contains 0 elements based on the search statement the fasta file will not be created
    def create_fasta(self):
        search_handle = Entrez.esearch(db="nucleotide",
                                   term = self.search_statement,
                                   usehistory="y",
                                   idtype="acc",
                                   retmax = 10000)
        search_record = Entrez.read(search_handle)


        fetch_handle = Entrez.efetch(db = 'nucleotide',
                                 id = search_record["IdList"],
                                 rettype='gb')
        check = len(search_record['IdList'])

        if check != 0:
            #print("I go through here)")
            os.system(f"mkdir /home/projects/animalia_mito/data/{self.genus}")

            fetch_record = SeqIO.parse(fetch_handle, "gb")

            fetch_record = list(fetch_record)

            SeqIO.write(fetch_record, self.path, "fasta")

        
        return check

    
    def faster_reader(self):

        with open(self.path) as handle:
            for record in SeqIO.FastaIO.FastaIterator(handle):
                yield record

    def sort_fasta(self):
                
        best_ref = {}
        for entry in self.faster_reader():
            nCount = entry.seq.lower().count('n')
            best_ref.update({entry.id : nCount})
    

        id_sorted = dict(sorted(best_ref.items(), key = lambda item: item[1]))

        end_sort = []
        for i in id_sorted.keys():
            for entry in self.faster_reader():
                if i == entry.id:
                    end_sort.append(entry)
        
            
        SeqIO.write(list(end_sort), self.path, "fasta" )
            
            
    def exclude_samples(self):
    
        kept = []
        seen = set()
        for entry in self.faster_reader():
            sp_name = entry.description.split()
            subsp_name = sp_name[1:4]
            if subsp_name[-1] == "mitochondrion,":
                sp = sp_name[1:3]
                sp = ' '.join(sp)
                if sp not in seen:
                    seen.add(sp)
                    kept.append(entry)
            else:
                subsp = ' '.join(subsp_name)
                if subsp not in seen:
                    seen.add(subsp)
                    kept.append(entry)
            
            
        SeqIO.write(list(kept), self.path2, "fasta" )
        
    def faster_reader2(self):

        with open(self.path2) as handle:
            for record in SeqIO.FastaIO.FastaIterator(handle):
                yield record


    def make_lineage_file(self):
        for entry in self.faster_reader2():
            os.system(f"""esearch -db nucleotide -query "{entry.id}" | esummary | xtract -pattern TaxId -element TaxId >> /home/projects/animalia_mito/data/{self.genus}/tax_id.txt""")


        with open(f'/home/projects/animalia_mito/data/{self.genus}/tax_id.txt') as index:
            test = index.read()

        tax_list = test.split('\n')
        tax_list = tax_list[:-1]


        lineage = {}
        id_list = []
        for entry in self.faster_reader2():
            id_list.append(entry.id)
        lin_list = []
        for i in tax_list:
            v = self.ncbi.get_lineage(i)
            names = self.ncbi.get_taxid_translator(v)
            a = [names[taxid] for taxid in v]
            lin_list.append(a)

        for i in range(len(id_list)):
            lineage[id_list[i]] = lin_list[i]

        output = open(f"/home/projects/animalia_mito/data/{self.genus}/{self.genus}_lin.txt", "w")
        for key, value in lineage.items():
            output.write('%s:%s\n' % (key, value))

        output.close()
        with open(f"/home/projects/animalia_mito/data/{self.genus}/{self.genus}_lin.txt", "r") as file:
            contents = file.read()
            n = 0
            if "Arthopoda" or "Tetrapoda" in contents:
                #import pdb; pdb.set_trace()
                n = n + 1
            else:
                n = 0
            return n

