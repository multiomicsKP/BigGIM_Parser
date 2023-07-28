import pandas as pd
input_dir = "/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/TCGA_mut/"

#Generate an ID mapping dictionary for gene symbol and NCBI gene ID
ID_map = pd.read_csv("/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/ID_convert/Table_ID_convert_9.19.2021.tsv", sep = '\t')
dic_ID_symbol2ncbi ={}
dic_ID_ncbi2symbol = {}
import math
ID_map.astype(str)
for i in range(0,ID_map.shape[0]):
    symbol = ID_map['Approved_symbol'][i]
    
    ncbi_id = ID_map['NCBI_Gene_ID'][i]
    ensemble_id = ID_map['Ensembl_gene_ID'][i]
    if math.isnan(ncbi_id) == False :
        ncbi_id = int(ncbi_id)
        if symbol not in dic_ID_symbol2ncbi:
            dic_ID_symbol2ncbi[symbol] = ncbi_id
        else:
            print(symbol)
for i in range(0,ID_map.shape[0]):
    alias_list = str(ID_map["Alias_symbols"][i]).split(',')
    ncbi_id = ID_map['NCBI_Gene_ID'][i]
    ensemble_id = ID_map['Ensembl_gene_ID'][i]
    if len(alias_list) > 0:
        for alias in alias_list:
            if alias not in dic_ID_symbol2ncbi:
                dic_ID_symbol2ncbi[alias] = ncbi_id
dic_ID_symbol2ncbi['MEN1'] = 4221
dic_ID_symbol2ncbi['TSC1']= 7248
dic_ID_symbol2ncbi['WHSC1'] = 7468
dic_ID_symbol2ncbi['MET'] = 4233
dic_ID_symbol2ncbi['HIST1H1E'] = 3008
dic_ID_symbol2ncbi['FLNA'] = 2316
dic_ID_symbol2ncbi["FAM46D"] = 169966
dic_ID_symbol2ncbi["HIST1H1C"] = 3006

print(dic_ID_symbol2ncbi['HCCS'])

input_dir1 = "/Users/guangrong/Documents/GitHub_project/fastqpi_BigGIM/KGs/PPI/"


Biogrid_data = pd.read_csv(input_dir1 + "BIOGRID_Human.csv", low_memory=False)
Biogrid_data = Biogrid_data
print(set(Biogrid_data['Experimental System']))

subject_id_list = Biogrid_data['Entrez Gene Interactor A']
object_id_list = Biogrid_data['Entrez Gene Interactor B']

subject_name_list = Biogrid_data["Official Symbol Interactor A"]
object_name_list = Biogrid_data["Official Symbol Interactor B"]

print(len(subject_id_list))

for i in range(0,len(subject_id_list)):
    if subject_name_list[i] not in dic_ID_symbol2ncbi:
        print(subject_name_list[i] + " not in sub dic")
    else:
        if dic_ID_symbol2ncbi[(subject_name_list[i])] != int(subject_id_list[i]):
            print(subject_id_list[i] + subject_name_list[i])
    
    if object_name_list[i] not in dic_ID_symbol2ncbi:
        print(object_name_list[i] + " not in obj dic")

    else:
        if dic_ID_symbol2ncbi[(object_name_list[i])] != int(object_id_list[i]):
            print(object_id_list[i] + object_name_list[i])
        

subject_category_list = []
object_category_list = []

predicate_list = []
for predicate in Biogrid_data['Experimental System']:
    if predicate == 'physical':
        predicate_list.append("physically_interacts_with")
        subject_category_list.append("Protein")
        object_category_list.append("Protein")
    elif predicate == "genetic":
        predicate_list.append("genetically_interacts_with")
        subject_category_list.append("Gene")
        object_category_list.append("Gene")

publications_list = Biogrid_data['Author']
publications_list_new = []
for p in publications_list:
    if p.startswith("PUBMED"):
        publications_list_new.append("PMID:" + str(p.split(":")[1]))
    else:
        publications_list_new.append(p)






Result = pd.DataFrame()
Result['subject_id'] = subject_id_list
Result['object_id'] = object_id_list
Result['subject_id_prefix'] = ['NCBIGene'] * len(subject_id_list)
Result['object_id_prefix'] = ['NCBIGene'] * len(object_id_list)
Result['subject_name'] = subject_name_list
Result['object_name'] = object_name_list
Result['predicate'] = predicate_list
Result['knowledge_source'] = ['Biogrid'] * len(subject_id_list)
Result['publications'] = publications_list_new
Result['subject_category'] = subject_category_list
Result['object_category'] = object_category_list

Result_filter = Result.loc[Result['object_name'].isin(dic_ID_symbol2ncbi.keys())]
Result_filter = Result_filter.loc[Result_filter['subject_name'].isin(dic_ID_symbol2ncbi.keys())]

Result_filter.to_csv("../Input_Table/Biogrid_formated_new.csv", index = None)