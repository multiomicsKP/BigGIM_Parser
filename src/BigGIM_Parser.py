# Author: Guangrong Qin 
# Institute for Systems Biology
# Aug 16, 2022

# reference of biolink model:
# https://github.com/biolink/biolink-model/blob/master/biolink-model.yaml
# infores CURIE: https://docs.google.com/spreadsheets/d/1Ak1hRqlTLr1qa-7O0s5bqeTHukj9gSLQML1-lg6xIHM/edit#gid=293462374    
# KG standarization: https://docs.google.com/document/d/1SlqN7M1LTgwBcoIuEAlUnxkJ3xzTmKJJ327zO9Vw5uQ/edit

import validators    
import pandas as pd
import numpy as np

#############Part I: KG overview ####################
def graph_stat(G_G_KG_formated):
    num_subject_ids = len(set(G_G_KG_formated['subject_id']))
    num_object_ids = len(set(G_G_KG_formated['object_id']))
    num_edges = G_G_KG_formated.shape[0]
    num_unique_nodes = len(set(list(G_G_KG_formated['subject_id']) + list(G_G_KG_formated['object_id'])))
    num_unique_association = len(set(G_G_KG_formated['predicate']))
    
    dic_node_degree = {}
    node_degree_list = []
    node_list = list(set(list(set(G_G_KG_formated['subject_id'])) + list(set(G_G_KG_formated['object_id']))))
    for node in node_list:
        x1 = list(G_G_KG_formated.loc[G_G_KG_formated['subject_id'] == node ].index)
        x2 = list(G_G_KG_formated.loc[G_G_KG_formated['object_id'] == node ].index)
        x=set(x1 + x2)
        node_degree_list.append(len(x))
        dic_node_degree[node] = len(x)
        
    node_degree_max = np.max(node_degree_list)
        
    print("max degree: " + str(node_degree_max))
    print("num of subject ids:" + str(num_subject_ids))
    print("num of object ids:" + str(num_object_ids))
    print("number of edges:" + str(num_edges))
    print("number of unique nodes: " + str(num_unique_nodes))
    print("number of association categories: " + str(num_unique_association))

#############Part II: KG parser ####################

def get_xref(id_prefix,id):
    
    dic_xref = {
            'NCBIGene':"https://www.ncbi.nlm.nih.gov/gene/",
            'Pubchem.compound': "https://pubchem.ncbi.nlm.nih.gov/compound/",
            'MONDO':"http://purl.obolibrary.org/obo/MONDO_",
           }

    if id_prefix in dic_xref:
        return(dic_xref[id_prefix]  + str(id))
    else:
        print("id_prefix not recognized!")
        return(None)


def header_check(file_formated):
    column_names = file_formated.columns.values.tolist()
    # minumus columns: 
    #'subject_id', 'subject_category', 'subject_name', 'subject.id prefixes', 
    #'object_id', 'object_category', 'object_name', 'object.id prefixes', 
    #'predicate', 'Knowledge_source']
    #['P_value']
    
    if "subject_id" in column_names and "subject_category" in column_names and "subject_name" in column_names and 'subject_id_prefixes' in column_names and 'object_id' in column_names and 'object_category' in column_names and 'object_name' and column_names and 'object_id_prefixes' in column_names and 'predicate' in column_names:
        format_checker = True
    else:
        format_checker = False
        print("Need the essential components: eg. subject_id, subject_category, subject_name, subject.id prefixes, object_id, object_category, object_name, object.id prefixes,predicate")
    return(format_checker)

    
def load_tsv_data(file_formated, Date, title):
    ## Pre-define associations
    correlation_statistic = {
        "T_test": {
                "attribute_type_id": "NCIT:C53236",  # Correlation Test -- http://purl.obolibrary.org/obo/NCIT_C53236
                "description": "t-test was used to compute the p-value for the association",
                "value": "NCIT:C53231",  # t-Test -- http://purl.obolibrary.org/obo/NCIT_C53231
                "value_type_id": "biolink:id"
        },
        "Spearman_correlation": {
                "attribute_type_id": "NCIT:C53236",  # Correlation Test -- http://purl.obolibrary.org/obo/NCIT_C53236
                "description": "Spearman Correlation Test was used to extract the association",
                "value": "NCIT:C53249",  # Spearman Correlation Test -- http://purl.obolibrary.org/obo/NCIT_C53249
                "value_type_id": "biolink:id",
        },
        "Wilcoxon-test":{
                "attribute_type_id": "NCIT:C53246",
                "description": "Wilcoxon-test was used to compute the p-value for the association",
                "value":"NCIT:C53246",
                "value_type_id": "biolink:id",
        },
        "Pearson_correlation":{
                "attribute_type_id": "NCIT:C53244", 
                "description": "Pearson correlation was used to extract the association",
                "value":"NCIT:C53246",  # Pearson correlation test -- https://ontobee.org/ontology/NCIT?iri=http://purl.obolibrary.org/obo/NCIT_C53244
                "value_type_id": "biolink:id",
        }}

    column_names = file_formated.columns.values.tolist()
    
    format_checker = header_check(file_formated)
        
        
    if format_checker == True:
        # ID validation
        validated_id = set()


        for index, row in file_formated.iterrows():
            unique_id_list = ["Multiomics-BigGIM-DrugResponse",str(row["subject_id"]),str(row["object_id"]), row['predicate'],row['knowledge_source']]
           
            #standerize id prefix
            subject_id_prefix = row['subject_id_prefixes']
            object_id_prefix = row['object_id_prefixes']
            subject_id_prefix_columnname = '_'.join(subject_id_prefix.split('.'))
            object_id_prefix_columnname = '_'.join(object_id_prefix.split('.'))
            
            # standarized ids for object and subject
            subject_id = subject_id_prefix + ":" + str(row["subject_id"])
            raw_subject_id = str(row["subject_id"])
            object_id = object_id_prefix + ":" + str(row["object_id"])
            raw_object_id = str(row["object_id"])
            
            # Get the dictionary of subjects, objects and associations
            
            subjects = {
                        "id": subject_id,
                        subject_id_prefix_columnname:subject_id,
                        "name": row["subject_name"],
                        "type": row["subject_category"], 
                        "xref": get_xref(row['subject_id_prefixes'],raw_subject_id)
                    }
            
            if subjects["xref"] not in validated_id:
                if validators.url(subjects["xref"]) == True:
                    validated_id.append(validators.url(subjects["xref"]))
                else:
                    print(subjects["xref"] + " is not valid")

            objects = {
                        "id": object_id,
                        object_id_prefix_columnname:object_id,
                        "name": row["object_name"],
                        "type": row["object_category"],
                        "xref": get_xref(row['object_id_prefixes'],raw_object_id)
                        }
            
            if objects["xref"] not in validated_id:
                if validators.url(objects["xref"]) == True:
                    validated_id.append(validators.url(objects["xref"]))
                else:
                    print(objects["xref"] + " is not valid. " )

            predicates = "biolink:"+'_'.join(row['predicate'].split(' '))
            
            
            edge_attributes = []
            
            association = {'edge_label':predicates,
                          "edge_attributes":edge_attributes
                          }
            
            
                
            # aggregator_knowledge_source
            edge_attributes.append({"attribute_type_id": "biolink:aggregator_knowledge_source",
                                    "value": "infores:biothings-multiomics-biggim-drugresponse"})
            
            #creation_date
           # edge_attributes.append({"attribute_type_id": "biolink:creation_date","value": Date})
            
            #resource_url
            edge_attributes.append({"attribute_type_id": "biolink:supporting_study_method_description", 
                                     "value": "https://github.com/NCATSTranslator/Translator-All/wiki/MultiomicsBigGIM_KP"})
            
            
            #P_value
            if "statistics_method" in column_names:
                attributes = correlation_statistic[row['statistics_method']]
                
            if "P_value" in column_names:
                edge_attributes.append(
                {
                    "attribute_type_id": "EDAM:data_0951",  # statistical estimate score -- http://edamontology.org/data_0951
                    #"description": "Confidence metric for the association",
                    "value": float(row["P_value"]),
                    "value_type_id": "EDAM:data_1669",  # P-value -- http://edamontology.org/data_1669
                    "attributes": [attributes]  # sub-attributes should be a list per TRAPI standard
                })
            
            
                
            #Sample size
            if "supporting_study_size" in column_names:
                edge_attributes.append(
                    {
                        "attribute_type_id": "biolink:supporting_study_size",  
                        #"description": "The sample size used in a study that provided evidence for the association",
                        "value": int(row['supporting_study_size']),
                    }
                )
                
            
            #Datasets for extracting the knowledge graphs
            if "Data_set" in column_names:
                dataset_attributes=[] # sub-attributes should be a list per TRAPI standard
                if 'PMID' in row["Data_set"]:
                    dataset_attributes.append(
                        {
                            "attribute_type_id": "biolink:Publication",
                          #  "description": "Publication describing the dataset used to compute the association",
                            "value": row["Data_set"],
                          #  "value_type_id": "biolink:id"
                        })
                    
                elif validators.url(row["Data_set"]) :
                     dataset_attributes.append(
                            {
                            "attribute_type_id": "biolink:source_url",
                          #  "description": "source url describing the dataset used to compute the association",
                            "value": row["Data_set"],
                          #  "value_type_id": "biolink:id"
                            })
                elif row["Data_set"].startswith("www.") and validators.url("https://"+row["Data_set"]):
                    dataset_attributes.append(
                            {
                            "attribute_type_id": "biolink:source_url",
                           # "description": "source url describing the dataset used to compute the association",
                            "value": "https://"+row["Data_set"],
                           # "value_type_id": "biolink:id"
                            })
                elif row["Data_set"].startswith("www.") and validators.url("http://"+row["Data_set"]):
                    dataset_attributes.append(
                            {
                            "attribute_type_id": "biolink:source_url",
                           # "description": "source url describing the dataset used to compute the association",
                            "value": "http://" +row["Data_set"],
                           # "value_type_id": "biolink:id"
                            })
                else:
                        dataset_attributes.append({"attribute_type_id":None,
                                                 "value":row['Data_set']})
                
                edge_attributes.append(
                {   "attribute_type_id": "biolink:Dataset",
                    #"description": "Dataset used to compute the association",
                     "value":row['Data_set'],
                    "dataset_attributes":dataset_attributes})
            
            # knowledge graphs for extract the association
            if "knowledge_source" in column_names:
                dataset_attributes=[] # sub-attributes should be a list per TRAPI standard
                if 'PMID' in row["knowledge_source"]:
                    dataset_attributes.append(
                        {
                            "attribute_type_id": "biolink:Publication",
                         #   "description": "Publication describing the association",
                            "value": row["knowledge_source"],
                         #   "value_type_id": "biolink:id"
                        }
                        )
                    
                elif validators.url(row["knowledge_source"]) :    
                    dataset_attributes.append(
                            {
                            "attribute_type_id": "biolink:source_url",
                          #  "description": "source url describing the association",
                            "value": row["knowledge_source"],
                          #  "value_type_id": "biolink:id"
                            })
                    
                elif row["knowledge_source"].startswith("www.") and validators.url("https://"+row["knowledge_source"]):
                    dataset_attributes.append(
                            {
                            "attribute_type_id": "biolink:source_url",
                           # "description": "source url describing the association",
                            "value": "https://"+row["knowledge_source"],
                           # "value_type_id": "biolink:id"
                            })
                elif row["knowledge_source"].startswith("www.") and validators.url("http://"+row["knowledge_source"]):
                    dataset_attributes.append(
                            {
                            "attribute_type_id": "biolink:source_url",
                           # "description": "source url describing association",
                            "value": "http://" +row["knowledge_source"],
                           # "value_type_id": "biolink:id"
                            })
                elif row["knowledge_source"] == 'DrugCentral':
                    dataset_attributes.append(
                        {
                            "attribute_type_id": "biolink:source_infores",
                          #  "description": "source infores describing association",
                            "value": "infores:drugcentral",
                          #  "value_type_id": "biolink:id"
                        })
                else:
                    dataset_attributes.append(
                        {
                            "attribute_type_id": None,
                            "value": row["knowledge_source"]
                        })
                        
                edge_attributes.append(
                {
                    "attribute_type_id": "biolink:knowledge_source",
                    #"value": row['knowledge_source'],
                    #"value_type_id": None,
                    "dataset_attributes":dataset_attributes
                })
                                 
            #add more optional associations
            # 

            #
            
            #Qualifiers
            if "subject_aspect_qualifier" in column_names:
                edge_attributes.append({
                    "attribute_type_id":"biolink:subject_aspect_qualifier",  # biolink version 3.0.0
                    "value":row["subject_aspect_qualifier"]
                })
                
            if "object_aspect_qualifier" in column_names:
                edge_attributes.append({
                    "attribute_type_id":"biolink:object_aspect_qualifier",  # biolink version 3.0.0
                    "value":row["object_aspect_qualifier"]
                })
                
             # disease context
            if "context_qualifier" in column_names:
                edge_attributes.append(
                {
                    "attribute_source": "infores:biothings-multiomics-biggim-drugresponse",
                    "attribute_type_id": "biolink:context_qualifier",  # biolink version 3.0.0
                    "value": row['context_qualifier'],
                    #"value_type_id": "biolink:id",
                })
            
            #anatomical context qualifier
            if "anatomical_context_qualifier" in column_names:
                edge_attributes.append(
                {
                    "attribute_type_id": "biolink:anatomical_context_qualifier",  # biolink version 3.0.0
                    "value": row['anatomical_context_qualifier'],
                })
                
            #frequence qualifier
            if "frequency_qualifier" in column_names:
                edge_attributes.append(
                {
                    "attribute_type_id": "biolink:frequency_qualifier",  # biolink version 3.0.0
                    "value": row['frequency_qualifier'],
                })
                
            json = {
                    #"_id":'-'.join(unique_id_list),
                    "_id":title + "_"+str(Date) + "_" + str(index),
                    "subject": subjects,
                    "association": association,
                    "object": objects
                    }
            yield json
