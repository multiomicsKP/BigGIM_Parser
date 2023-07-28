# Author: Guangrong Qin
# Institute for Systems Biology
# Aug 16, 2022

# reference of biolink model:
# https://github.com/biolink/biolink-model/blob/master/biolink-model.yaml
# infores CURIE: https://docs.google.com/spreadsheets/d/1Ak1hRqlTLr1qa-7O0s5bqeTHukj9gSLQML1-lg6xIHM/edit#gid=293462374
# KG standarization: https://docs.google.com/document/d/1SlqN7M1LTgwBcoIuEAlUnxkJ3xzTmKJJ327zO9Vw5uQ/edit

import os
import validators
import pandas as pd


############# Part I: KG overview ####################
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
        x1 = list(G_G_KG_formated.loc[G_G_KG_formated['subject_id'] == node].index)
        x2 = list(G_G_KG_formated.loc[G_G_KG_formated['object_id'] == node].index)
        x = set(x1 + x2)
        node_degree_list.append(len(x))
        dic_node_degree[node] = len(x)

    node_degree_max = max(node_degree_list)

    print("max degree: " + str(node_degree_max))
    print("num of subject ids:" + str(num_subject_ids))
    print("num of object ids:" + str(num_object_ids))
    print("number of edges:" + str(num_edges))
    print("number of unique nodes: " + str(num_unique_nodes))
    print("number of association categories: " + str(num_unique_association))


############# Part II: KG parser ####################

# Pre-define associations
CORRELATION_STATISTIC = {
    "T_test": {
        # Correlation Test -- http://purl.obolibrary.org/obo/NCIT_C53236
        "attribute_type_id": "NCIT:C53236",
        "description": "t-test was used to compute the p-value for the association",
        "value": "NCIT:C53231",  # t-Test -- http://purl.obolibrary.org/obo/NCIT_C53231
        "value_type_id": "biolink:id"
    },
    "Spearman_correlation": {
        # Correlation Test -- http://purl.obolibrary.org/obo/NCIT_C53236
        "attribute_type_id": "NCIT:C53236",
        "description": "Spearman Correlation Test was used to extract the association",
        # Spearman Correlation Test -- http://purl.obolibrary.org/obo/NCIT_C53249
        "value": "NCIT:C53249",
        "value_type_id": "biolink:id",
    },
    "Wilcoxon-test": {
        "attribute_type_id": "NCIT:C53246",
        "description": "Wilcoxon-test was used to compute the p-value for the association",
        "value": "NCIT:C53246",
        "value_type_id": "biolink:id",
    },
    "Pearson_correlation": {
        "attribute_type_id": "NCIT:C53244",
        "description": "Pearson correlation was used to extract the association",
        # Pearson correlation test -- https://ontobee.org/ontology/NCIT?iri=http://purl.obolibrary.org/obo/NCIT_C53244
        "value": "NCIT:C53246",
        "value_type_id": "biolink:id",
    }
}


DIC_XREF = {
    'NCBIGene': "https://www.ncbi.nlm.nih.gov/gene/",
    'Pubchem.compound': "https://pubchem.ncbi.nlm.nih.gov/compound/",
    'PUBCHEM.COMPOUND': "https://pubchem.ncbi.nlm.nih.gov/compound/",
    'pubchem.compound': "https://pubchem.ncbi.nlm.nih.gov/compound/",
    'pubchem.cid': "https://pubchem.ncbi.nlm.nih.gov/compound/",
    'MONDO': "http://purl.obolibrary.org/obo/MONDO_",
    "CHEBI": "https://www.ebi.ac.uk/chebi/chebiOntology.do?chebiId=CHEBI:",
    "CHEMBL": "https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL",
    "ENSEMBL": "https://uswest.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",
    "NCIT": "https://ontobee.org/ontology/NCIT?iri=http://purl.obolibrary.org/obo/NCIT_",
    "CellOntology":"http://purl.obolibrary.org/obo/"
}


def get_xref(id_prefix, _id):
    if id_prefix in DIC_XREF:
        if id_prefix == "MONDO":
            return DIC_XREF[id_prefix] + str(_id).split(":")[1]  #MONDO ID in the table is labeled as MONDO:0019391
        else:
            return DIC_XREF[id_prefix] + str(_id)
            
    else:
        print("id_prefix not recognized!")
        return None

def validate_xref(validated_xrefs, xref):
    if xref not in validated_xrefs:
        if validators.url(xref):
            validated_xrefs.add(xref)
        else:
            print(f"{xref} is not valid. ")


def header_check(file_formated):
    file_columns = set(file_formated.columns.values.tolist())

    # minumus columns:
    # 'subject_id', 'subject_category', 'subject_name', 'subject_id_prefix',
    # 'object_id', 'object_category', 'object_name', 'object_id_prefix',
    # 'predicate', 'Knowledge_source', 'P_value'
    essential_columns = [
        "subject_id",
        "subject_category",
        "subject_name",
        "subject_id_prefix",

        "object_id",
        "object_category",
        "object_name",
        "object_id_prefix",

        "predicate",
        # "knowledge_source"
    ]

    for col in essential_columns:
        if col not in file_columns:
            print(f"Missing essential column {col}")
            return False

    return True


def _parse_party(row, party):
    """For given input row and party, generate the subject/object json 

    Args:
        row: a row of a dataframe
        party: either "object" or "subject" (no error handling)
    """
    prefix_column = f"{party}_id_prefix"
    id_column = f"{party}_id"
    name_column = f"{party}_name"
    type_column = f"{party}_category"

    prefix = row[prefix_column]
    prefix_field_name = '_'.join(prefix.split('.'))

    raw_id = row[id_column]
    if type(raw_id) == float:
        raw_id = int(raw_id)
    raw_id = str(raw_id)

    if prefix == "MONDO" or prefix == "CellOntology":
        _id = raw_id
    else:
        _id = f"{prefix}:{raw_id}"

    xref = get_xref(prefix, raw_id)

    json = {
        "id": _id,
        prefix_field_name: _id,
        "name": row[name_column],
        "type": row[type_column],
        "xref": xref
    }
    return json


def parse_subject(row):
    return _parse_party(row, "subject")


def parse_object(row):
    return _parse_party(row, "object")


def parse_sub_attribute(source, infores_dict):
    """Parse sub-attributes for column "knowledge_source" or "Data_set"

    Args:
        source (str): the content of column "knowledge_source" or "Data_set"
        infores_dict (dict): a dict of infores keys (to match against source) and values (as part of the returned attributes). E.g.
                             {"GTEx": "infores:gtex", "HuRI": "infores:HuRI"}
    """
    def _parse_publication(source):
        if 'PMID' in source:
            new_values = source.split(":")
            new_value_PMID = new_values[0].strip() + ":" + new_values[1].strip()

            attribute = {
                # "description": "Publication describing the dataset used to compute the association",
                # "value_type_id": "biolink:id"
                "attribute_type_id": "biolink:Publication",
                "value": new_value_PMID,
            }
            return attribute

        return None

    def _parse_source_url(source):
        url = None
        if validators.url(source):
            url = source
        elif source.startswith("www."):
            for candidate_url in [f"https://{source}", f"http://{source}"]:
                if validators.url(candidate_url):
                    url = candidate_url
                    break

        if url is not None:
            attribute = {
                # "description": "source url describing the dataset used to compute the association",
                # "value_type_id": "biolink:id"
                "attribute_type_id": "biolink:source_url",
                "value": url
            }
            return attribute

        return None

    def _parse_source_infores(source, infores_dict):
        for infores_key, infores_value in infores_dict.items():
            if source == infores_key:
                attribute = {
                    # "description": "source infores describing association",
                    # "value_type_id": "biolink:id"
                    "attribute_type_id": "biolink:source_infores",
                    "value": infores_value,
                }
                return attribute

        return None

    def _default_attribute(source):
        return {"attribute_type_id": "Biolink:Dataset", "value": source}

    # main function body starts
    # return the first non-None value
    attribute = _parse_publication(source) or _parse_source_url(source) or \
        _parse_source_infores(source, infores_dict) or _default_attribute(source)
    return attribute


def parse_edge_attributes(row, column_names):
    edge_attributes = []

    # aggregator_knowledge_source
    # edge_attributes.append({
    #    "attribute_type_id": "biolink:aggregator_knowledge_source",
    #    "value": "infores:biothings-multiomics-biggim-drugresponse"
    # })

    # creation_date
    # edge_attributes.append({"attribute_type_id": "biolink:creation_date","value": Date})

    # resource_url
    # edge_attributes.append({"attribute_type_id": "biolink:supporting_study_method_description",
    #                         "value": "https://github.com/NCATSTranslator/Translator-All/wiki/MultiomicsBigGIM_KP"})

    # P_value
    if "statistics_method" in column_names and "P_value" in column_names:
        edge_attributes.append({
            # statistical estimate score -- http://edamontology.org/data_0951
            "attribute_type_id": "EDAM:data_0951",
            # "description": "Confidence metric for the association",
            "value": float(row["P_value"]),
            "value_type_id": "EDAM:data_1669",  # P-value -- http://edamontology.org/data_1669
            # sub-attributes should be a list per TRAPI standard
            "attributes": [CORRELATION_STATISTIC[row['statistics_method']]]
        })

    # Sample size
    if "supporting_study_size" in column_names:
        edge_attributes.append({
            "attribute_type_id": "biolink:supporting_study_size",
            # "description": "The sample size used in a study that provided evidence for the association",
            "value": int(row['supporting_study_size']),
        })

    # Datasets for extracting the knowledge graphs
    if "Data_set" in column_names:
        source = row["Data_set"]
        # infores_dict = {"GTEx": "infores:gtex"}
        # attribute = parse_sub_attribute(source, infores_dict)

        edge_attributes.append({
            # "description": "Dataset used to compute the association",
            "attribute_type_id": "biolink:dataset",
            # "attribute_type_id": "biolink:Data_source",
            "value": source,
            # "attributes": [attribute]  # sub-attributes should be a list per TRAPI standard
        })

    # publications
    if "publications" in column_names:
        publications = row["publications"]
        if 'PMID' in publications:
            new_values = publications.split(":")
            new_value_PMID = new_values[0].strip() + ":" + new_values[1].strip()

            edge_attributes.append({
                "attribute_type_id": "biolink:publications",
                "value": new_value_PMID,
            })
        else:
            edge_attributes.append({
                "attribute_type_id": "biolink:publications",
                "value": publications,
            })

    # knowledge graphs for extract the association
    if "knowledge_source" in column_names:
        source = row["knowledge_source"]
        infores_dict = {
            "Biogrid": "infores:biogrid",
            "HuRI": "infores:huri",
            "DrugCentral": "infores:drugcentral",
            "http://www.interactome-atlas.org/download": "infores:huri",
            "TTD_2021": "infores:ttd",
            "CellMarker": "infores:cellmarker2.0",
        }

        attribute = parse_sub_attribute(source, infores_dict)
        if source in infores_dict:
            edge_attributes.append({
                "attribute_type_id": "biolink:primary_knowledge_source",
                # "value": source,
                "value": infores_dict[source], 
                # "value_type_id": None,
                "attributes": [attribute]  # sub-attributes should be a list per TRAPI standard
            })

            edge_attributes.append({
                "attribute_type_id": "biolink:aggregator_knowledge_source",
                "value": "infores:biothings-multiomics-biggim-drugresponse"
            })

        elif source == "BigGIM_DRUGRESPONSE":
            edge_attributes.append({
                "attribute_type_id": "biolink:primary_knowledge_source",
                "value": "infores:biothings-multiomics-biggim-drugresponse"
            })

        elif source == "BigGIM":
            edge_attributes.append({
                "attribute_type_id": "biolink:primary_knowledge_source",
                "value": "infores:biggim"
            })
        elif source == "CellMarker":
            edge_attributes.append({
                "attribute_type_id": "biolink:primary_knowledge_source",
                "value": "infores:cellmarker2.0"
            })
        elif "PMID" in source:
            edge_attributes.append({
                "attribute_type_id": "biolink:primary_knowledge_source",
                "value": "infores:pubmed"
            })
            edge_attributes.append({
                "attribute_type_id": "biolink:publications",
                "value": source
            })

    # add more optional associations

    # Qualifiers
    if "subject_aspect_qualifier" in column_names:
        edge_attributes.append({
            "attribute_type_id": "biolink:subject_aspect_qualifier",  # biolink version 3.0.0
            "value": row["subject_aspect_qualifier"]
        })

    if "object_aspect_qualifier" in column_names:
        edge_attributes.append({
            "attribute_type_id": "biolink:object_aspect_qualifier",  # biolink version 3.0.0
            "value": row["object_aspect_qualifier"]
        })

    # disease context
    if "context_qualifier" in column_names:
        edge_attributes.append({
            "attribute_source": "infores:biothings-multiomics-biggim-drugresponse",
            "attribute_type_id": "biolink:context_qualifier",  # biolink version 3.0.0
            "value": row['context_qualifier'],
            # "value_type_id": "biolink:id",
        })

    # anatomical context qualifier
    if "anatomical_context_qualifier" in column_names:
        edge_attributes.append({
            "attribute_type_id": "biolink:anatomical_context_qualifier",  # biolink version 3.0.0
            "value": row['anatomical_context_qualifier'],
        })

    # frequence qualifier
    if "frequency_qualifier" in column_names:
        edge_attributes.append({
            "attribute_type_id": "biolink:frequency_qualifier",  # biolink version 3.0.0
            "value": row['frequency_qualifier'],
        })

    # supporting_study_cohort
    if "supporting_study_cohort" in column_names:
        edge_attributes.append({
            "attribute_type_id": "biolink:supporting_study_cohort",  # biolink version 3.0.0
            "value": row['supporting_study_cohort'],
        })

    # has_count
    if "has_count" in column_names:
        edge_attributes.append({
            "attribute_type_id": "biolink:has_count",  # biolink version 3.0.0
            "value": row['has_count'],
        })

    return edge_attributes


def parse_source_attribute(row, column_names):
    source_attributes = []

    infores_dict = {
        "Biogrid": "infores:biogrid",
        "HuRI": "infores:huri",
        "DrugCentral": "infores:drugcentral",
        "http://www.interactome-atlas.org/download": "infores:huri",
        "TTD_2021": "infores:ttd",
        "TCGA": "infores:tcga",
        "GDSC": "infores:gdsc",
        "GTEx": "infores:gtex",
        "CellMarker": "infores:cellmarker2.0"
    }
    
    if "knowledge_source" in column_names:
        source = row["knowledge_source"]

        # attribute = parse_sub_attribute(source, infores_dict)
        if source in infores_dict:
            source_attributes.append({
                "resource_id": infores_dict[source],
                "resource_role": "primary_knowledge_source"
            })

            source_attributes.append({
                "resource_id": "infores:biothings-multiomics-biggim-drugresponse",
                "resource_role": "aggregator_knowledge_source",
                "upstream_resource_ids": [infores_dict[source]]

            })

        elif source == "BigGIM_DRUGRESPONSE":
            source_attributes.append({
                "resource_id": "biolink:primary_knowledge_source",
                "resource_role": "primary_knowledge_source"
            })

        elif source == "BigGIM":
            source_attributes.append({
                "resource_id": "infores:biothings-multiomics-biggim-drugresponse",
                "resource_role": "primary_knowledge_source"
            })
        
        elif "PMID" in source:
            source_attributes.append({
                "resource_id": "infores:pubmed",
                "resource_role": "primary_knowledge_source"
            })
            source_attributes.append({
                "resource_id":  "infores:biothings-multiomics-biggim-drugresponse",
                "resource_role": "aggregator_knowledge_source",
                "upstream_resource_ids": source
            })
        else:
            source_attributes.append({
                "resource_id": "infores:biothings-multiomics-biggim-drugresponse",
                "resource_role": "primary_knowledge_source"
            })
    else:
        source_attributes.append({
                "resource_id": "infores:biothings-multiomics-biggim-drugresponse",
                "resource_role": "primary_knowledge_source"
            })
        if "Data_set" in column_names:
            source = row["Data_set"]
            # attribute = parse_sub_attribute(source, infores_dict)

            if source in infores_dict:
                source_attributes.append({
                    "resource_id": infores_dict[source],
                    "resource_role": "supporting_data_source"
                })

            elif "PMID" in source:
                source_attributes.append({
                    "resource_id": "infores:pubmed",
                    "resource_role": "supporting_data_source",
                  #  "resource_name": source
                })

        if "supporting_study_cohort" in column_names:
            if row['supporting_study_cohort'] in infores_dict:
                source_attributes.append({
                    "resource_id": infores_dict[row['supporting_study_cohort']],
                    "resource_role": "supporting_data_set"} ) # biolink version 3.0.0
            else:
                print("supporting_study_cohort not in infores_dict")

    return source_attributes


def load_tsv_data(file_path):
    file_formated = pd.read_csv(file_path)
    file_formated = file_formated.dropna()
    if not header_check(file_formated):
        print(f"file {file_path} misses essential headers.")
        return

    column_names = file_formated.columns.values.tolist()
    validated_xrefs = set()

    for index, row in file_formated.iterrows():
        # Get the dictionary of subject, object and associations

        subject_json = parse_subject(row)
        validate_xref(validated_xrefs, subject_json["xref"])

        object_json = parse_object(row)
        validate_xref(validated_xrefs, object_json["xref"])

        predicates = "biolink:" + '_'.join(row['predicate'].split(' '))
        edge_attributes = parse_edge_attributes(row, column_names)
        source_attributes = parse_source_attribute(row, column_names)

        association_json = {
            "edge_label": predicates,
            "edge_attributes": edge_attributes,
            "source_attributes": source_attributes
        }

        json = {
            # "_id":'-'.join(unique_id_list),
            "_id": subject_json["type"] + "_" + predicates + "_" + object_json["type"] + "_" + file_path.split("/")[-1] + "_" + str(index),
            "subject": subject_json,
            "object": object_json,
            "predicate": association_json["edge_label"],
            "attributes": association_json["edge_attributes"],
            "sources": association_json["source_attributes"]  
        
        }
        yield json


def load_data(data_folder):
    file_names = [
        "Input_DrugResponse_expr_auc_gdsc_08312022.csv",
        "Input_DrugResponse_mut_IC50_gdsc_08312022.csv",
        "TCGA_driver_mut_freq.csv",
        "FA_mut.csv",
        "GDSC_cancer_specific_signatures.csv",
        "Drug_target_with_primary_source.csv",
        "Biogrid_formated.csv",
        "H-I-05_formated.csv",
        "HI-II-14_formated.csv",
        "HuRI_formated.csv",
        "Yang-16_formated.csv",
        "GTEX_liver_negative_correlated_formated.csv",
        "GTEX_liver_positively_correlated_formated.csv"
    ]
    file_paths = [os.path.join(data_folder, fn) for fn in file_names]
    for file_path in file_paths:
        for json in load_tsv_data(file_path):
            yield json
