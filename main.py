#Author: Guangrong_Qin
#Date: 2022.08.24

import json
import sys
sys.path.append("/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/src/")
from BigGIM_Parser import * 

#Datafile = "/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/formated/gene_gene_formatted.csv"
#TCGA_driver_mutations
Datafile1 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/TCGA_driver_mut_freq.csv"

Datafile2 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/FA_mut.csv"

def main(verbose = False):
    counter = 0
    for row in load_tsv_data(Datafile1,'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1
    #print(counter)
    for row in load_tsv_data(Datafile2,'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1

if __name__ == "__main__":
    main(verbose=True)