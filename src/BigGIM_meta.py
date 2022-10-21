import sys
import pandas as pd
import os

def Meta(KG):
    No_nodes = 0
    No_edges = 0
    No_record = 0 
    nodes = set()
    predicates = set()
    No_categories_nodes = set()

    result = {}

    data  = pd.read_csv(KG)
    subjects = set(data['subject_id'])
    objects = set(data['object_id'])
    predicates = set(data['predicate'])
    subject_id_prefix = set(data['subject_id_prefix'])
    object_id_prefix = set(data['object_id_prefix'])
    No_record = data.shape[0]
    No_objects = len(objects)
    No_subjects = len(subjects)

    No_nodes = len(set(list(subjects)+list(objects)))
    categories = set(list(data['subject_category']) + list(data['object_category']))
    result = {
            "KG":KG,
            "No_record":No_record,
            "No_nodes":No_nodes,
            "No_subject":No_subjects,
            "No_objects":No_objects,
            "No_unique_edges":len(predicates),
            "unique_edges":predicates,
            "No_categories_nodes": categories,
            "subject_id_prefix": subject_id_prefix,
            "object_id_prefix": object_id_prefix
            }

    return(result)


def Load_meta_list(data_folder):
    meta=[]
    filename_paths = [  os.path.join(data_folder,"Input_DrugResponse_expr_auc_gdsc_08312022.csv"),
                        os.path.join(data_folder,"Input_DrugResponse_mut_IC50_gdsc_08312022.csv"),
                        os.path.join(data_folder,"TCGA_driver_mut_freq.csv"),
                        os.path.join(data_folder,"FA_mut.csv"),
                        os.path.join(data_folder,"GDSC_cancer_specific_signatures.csv"),
                        os.path.join(data_folder,"Drug_targets_14806.csv"),
                        os.path.join(data_folder,"Biogrid_formated.csv"),
                        os.path.join(data_folder,"H-I-05_formated.csv"),
                        os.path.join(data_folder,"HI-II-14_formated.csv"),
                        os.path.join(data_folder,"HuRI_formated.csv"),
                        os.path.join(data_folder,"Yang-16_formated.csv"),
                        os.path.join(data_folder,"GTEX_liver_negative_correlated_formated.csv"),
                        os.path.join(data_folder,"GTEX_liver_positively_correlated_formated.csv"),
                    ]
    for filename_path in filename_paths:
        print(filename_path)
        result = Meta(filename_path)
        meta.append(result)
    return(meta)

def get_meta_df(meta):
    KG_list = []
    No_record = []
    No_nodes = []
    No_unique_edges = []
    unique_edges = []
    No_categories_nodes = []
    subject_id_prefix = []
    object_id_prefix  = []
    No_subjects = []
    No_objects = []

    for record in meta:
        KG_list.append(record['KG'])
        No_record.append(record['No_record'])
        No_nodes.append(record['No_nodes'])
        No_subjects.append(record['No_subject'])
        No_objects.append(record['No_objects'])
        No_unique_edges.append(record['No_unique_edges'])
        unique_edges.append(record['unique_edges'])
        No_categories_nodes.append(record['No_categories_nodes'])
        subject_id_prefix.append(record['subject_id_prefix'])
        object_id_prefix.append(record['object_id_prefix'])
    
    result =  pd.DataFrame({
        "KG":KG_list,
        "No_record": No_record,
        "No_nodes": No_nodes,
        "No_objects":No_objects,
        "No_subjects":No_subjects,
        "No_unique_edges": No_unique_edges,
        "unique_edges": unique_edges,
        "No_categories_nodes":No_categories_nodes,
        "subject_id_prefix":subject_id_prefix,
        "object_id_prefix":object_id_prefix
    })
    return(result)



def main():
    meta = Load_meta_list("/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/")
    result = (get_meta_df(meta))
    result.to_csv("/Users/guangrong/Documents/GitHub_project/BigGIM_Parser/Input_Table/meta.csv")

if __name__ == "__main__":
    main()