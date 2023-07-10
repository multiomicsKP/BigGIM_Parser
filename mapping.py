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
                "pubchem_compound": {
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
                "pubchem_compound": {
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
                "upstream_resource_ids": {
                    "normalizer": "keyword_lowercase_normalizer",
                    "type": "keyword"
                }
            }
        }
    }

    return mapping
