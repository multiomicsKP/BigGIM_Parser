#Author: Guangrong_Qin
#Date: 2022.08.24

import json
import sys
sys.path.append("/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/src/")
from BigGIM_Parser import * 

def main(verbose = False):
    counter = 0
    for row in load_tsv_data("/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/formated/gene_gene_formatted.csv",'08/16/2022'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1
    #print(counter)

if __name__ == "__main__":
    main(verbose=True)