#Author: Guangrong_Qin
#Date: 2022.08.24

import json
import sys
sys.path.append("/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/src/")
from BigGIM_Parser import * 

#Datafile = "/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/formated/gene_gene_formatted.csv"
#TCGA_driver_mutations
verbose1 = False
Datafile1 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/TCGA_driver_mut_freq.csv" #gene to disease

verbose2 = False
Datafile2 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/FA_mut.csv" #gene to disease

verbose3 = False
Datafile3 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/GDSC_cancer_specific_signatures.csv" # disease to gene

verbose4 = False
Datafile4 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/Biogrid_formated.csv" # protein to protein / gene to gene

verbose5 = False
Datafile5 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/HuRI_formated.csv" #protein to protein

verbose6 = False
Datafile6 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/HI-II-14_formated.csv" #protein to protein

verbose7 = False
Datafile7 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/H-I-05_formated.csv" #protein to protein

verbose8 = False
Datafile8 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/Yang-16_formated.csv" #protein to protein

verbose9 = False
Datafile9 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/GTEX_liver_negative_correlated_formated.csv" #gene expression negatively correlated

verbose10 = True
Datafile10 = "/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/GTEX_liver_positively_correlated_formated.csv" #gene expression negatively correlated

def main(verbose = False):
    counter = 0
    if verbose1:
        for row in load_tsv_data(Datafile1,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1

    if verbose2:
        for row in load_tsv_data(Datafile2,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1

    if verbose3:
        for row in load_tsv_data(Datafile3,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1
    
    if verbose4:
        for row in load_tsv_data(Datafile4,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1
    
    if verbose5:
        for row in load_tsv_data(Datafile5,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1

    if verbose6:
        for row in load_tsv_data(Datafile6,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1
    
    if verbose7:
        for row in load_tsv_data(Datafile7,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1
    
    if verbose8:
        for row in load_tsv_data(Datafile8,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1
    
    if verbose9:
        for row in load_tsv_data(Datafile9,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1
    
    if verbose10:
        for row in load_tsv_data(Datafile10,'08/16/2022'):
            print(json.dumps(row, sort_keys=True, indent=2))
            counter += 1

if __name__ == "__main__":
    main(verbose=True)