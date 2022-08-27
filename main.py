#Author: Guangrong_Qin
#Date: 2022.08.24

import json
import sys
sys.path.append("/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/src/")
from BigGIM_Parser import * 

#Datafile = "/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/formated/gene_gene_formatted.csv"
#TCGA_driver_mutations
Datafile1 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/TCGA_driver_mut_freq.csv" #gene to disease
  
Datafile2 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/FA_mut.csv" #gene to disease

Datafile3 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/GDSC_cancer_specific_signatures.csv" # disease to gene

Datafile4 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/Biogrid_formated.csv" # protein to protein / gene to gene

Datafile5 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/HuRI_formated.csv" #protein to protein

Datafile6 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/HI-II-14_formated.csv" #protein to protein

Datafile7 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/H-I-05_formated.csv" #protein to protein

def main(verbose = False):
    counter = 0
    for row in load_tsv_data(Datafile1,'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1
    for row in load_tsv_data(Datafile2,'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1

    for row in load_tsv_data(Datafile3,'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1
    
    #for row in load_tsv_data(Datafile4,'08/16/2022'):
    #    if verbose:
    #        print(json.dumps(row, sort_keys=True, indent=2))
    #    counter += 1

    for row in load_tsv_data(Datafile5,'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1

    for row in load_tsv_data(Datafile6,'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1
if __name__ == "__main__":
    main(verbose=True)