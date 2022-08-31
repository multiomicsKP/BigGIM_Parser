#Author: Guangrong_Qin
#Date: 2022.08.31

import json
import sys
sys.path.append("./src/") 
from BigGIM_Parser import * 

def main():
    counter = 0
    config_file = open("./Input_Table/config.file.txt")
    for line in config_file.readlines():
        #print(line)
        if line.startswith("#") == False and len(line)>0:
            filename   = line.split(',')[0].strip()
            date_label = line.split(',')[1].strip()
            label      = line.split(',')[2].strip()

            if label == "ON":
                for row in load_tsv_data(filename,date_label):
                    print(json.dumps(row, sort_keys=True, indent=2))
                    counter += 1
    
    config_file.close()

if __name__ == "__main__":
    main()