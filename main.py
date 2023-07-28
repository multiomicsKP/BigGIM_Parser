#Author: Guangrong_Qin
#Date: 2022.08.31

import json
import BigGIM_Parser 
from BigGIM_Parser import load_tsv_data
import os

#def main():
#    counter = 0
#    config_file = open("./Input_Table/config.file.txt")
#    for line in config_file.readlines():
#        #print(line)
#        if line.startswith("#") == False and len(line)>0:
#            filename   = line.split(',')[0].strip()
#            date_label = line.split(',')[1].strip()
#            label      = line.split(',')[2].strip()#
#
#            if label == "ON":
#                for row in load_tsv_data(filename,date_label):
#                    print(json.dumps(row, sort_keys=True, indent=2))
#                    counter += 1
#    
#    config_file.close()


def load_data(data_folder):
    filename_paths = [  #os.path.join(data_folder,"Input_DrugResponse_expr_auc_gdsc_08312022.csv"), #tested
                        #os.path.join(data_folder,"Input_DrugResponse_mut_IC50_gdsc_08312022.csv"), #tested
                        #os.path.join(data_folder,"TCGA_driver_mut_freq.csv"), #tested
                        #os.path.join(data_folder,"FA_mut.csv"), #tested
                        #os.path.join(data_folder,"GDSC_cancer_specific_signatures.csv"), #tested
                        #os.path.join(data_folder,"Drug_target_with_primary_source.csv"),
                        #os.path.join(data_folder,"Biogrid_formated.csv"), #tested
                        #os.path.join(data_folder,"H-I-05_formated.csv"), #tested
                        #os.path.join(data_folder,"HI-II-14_formated.csv"), #tested
                        #os.path.join(data_folder,"HuRI_formated.csv"), #tested
                        #os.path.join(data_folder,"Yang-16_formated.csv"),#tested
                        #os.path.join(data_folder,"GTEX_liver_negative_correlated_formated.csv"),#tested
                        #os.path.join(data_folder,"GTEX_liver_positively_correlated_formated.csv"), #tested
                        os.path.join(data_folder,"cellmarker.csv"), 
                    ]
    for filename_path in filename_paths:                
        for row in load_tsv_data(filename_path):
            yield row


def main():
    counter = 0
    verbose = True
    for row in load_data('Input_Table/'):
        if verbose:
            print(json.dumps(row, sort_keys=True, indent=2))
        counter += 1



if __name__ == "__main__":
    main()


