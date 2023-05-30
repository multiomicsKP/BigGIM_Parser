def custom_mapping(cls):
    """
    Method to provide custom mapping to parser.

    Configuration of this method in manifest.json should be:

        "uploader" : {
            "mapping" : "mapping:custom_mapping"
        }

    This is a class method but @classmethod decorator is not necessary. 
    See https://docs.biothings.io/en/latest/tutorial/studio_guide.html#manifest-based-data-plugins
    """
    mapping = {
        "subject": {
            "properties": {
                "id": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "name": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "type": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "xref": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "Pubchem_compound": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "MONDO": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "NCBIGene": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                }
            }
        },

         "predicate": {
            "normalizer": "keyword_lowercase_normalizer",
            "type": "keyword"
            },

        "attributes": {
            "properties": {
                "attribute_type_id": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "value_type_id": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "attributes": {
                    "properties": {
                        "attribute_type_id": {
                            "normalizer": "keyword_lowercase_normalizer",
                            "type": "keyword"
                        },
                        "value": {
                            "type": "text",
                        },
                        "value_type_id": {
                            "normalizer": "keyword_lowercase_normalizer",
                            "type": "keyword"
                        },
                        "description": {
                            "type": "text"
                        }
                    }
                },
                "attribute_source": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "value": {
                    "type": "text",
                    "index": False
                }
            }
        },

        "sources": {
            "properties": {
                "resource_id": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "resource_role": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                #"resource_name": {
                #    "normalizer": "keyword_lowercase_normalizer",
                #    "type": "keyword"
                #},
                "upstream_resource_ids": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                }
            }
        },
        
       # "association": {
       #     "properties": {
       #         "edge_label": {
       #             "normalizer": "keyword_lowercase_normalizer",
       #             "type": "keyword"
       #         },
       #         "edge_attributes": {
       #             "properties": {
       #                 "attribute_type_id": {
       #                     "normalizer": "keyword_lowercase_normalizer",
       #                     "type": "keyword"
       #                 },
       #                 "value_type_id": {
       #                     "normalizer": "keyword_lowercase_normalizer",
       #                     "type": "keyword"
       #                 },
       #                 "attributes": {
       #                     "properties": {
       #                         "attribute_type_id": {
       #                             "normalizer": "keyword_lowercase_normalizer",
       #                             "type": "keyword"
       #                         },
       #                         "value": {
       #                             "type": "text",
       #                             "index": False
       #                         },
       #                         "value_type_id": {
       #                             "normalizer": "keyword_lowercase_normalizer",
       #                             "type": "keyword"
       #                         },
       #                         "description": {
       #                             "type": "text"
       #                         }
       #                     }
       #                 },
       #                 "attribute_source": {
       #                     "normalizer": "keyword_lowercase_normalizer",
       #                     "type": "keyword"
       #                 },
       #                 "value": {
       #                     "type": "text",
       #                     "index": False
       #                 }
       #             }
       #         }
       #     }
       # },
        "object": {
            "properties": {
                "id": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "name": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "type": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "xref": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "NCBIGene": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "MONDO": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                },
                "PUBCHEM_COMPOUND": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                }
            }
        }
    }

    return mapping