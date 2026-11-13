from rnacentral_pipeline.databases.data import DiseaseAnnotation, Entry
from rnacentral_pipeline.databases.ols import search


def primary_id(record):
    return "EVLNCRNAS:" + record["ID"]


def accession(record):
    return "EVLNCRNAS:" + record["ID"]


def taxid(record):
    return int(record["taxid"])


def sequence(record):
    return record["sequence"].upper().replace("U", "T")


def aliases(record):
    if record["Aliases"] is None:
        return [str(record["external_id"])]

    aliases = [str(record["external_id"])]
    aliases.extend(a.strip() for a in str(record["Aliases"]).split(",") if a.strip())

    return aliases


def region_builder(info):
    return []


def rna_type(record):
    return record["so_type"]


def url(record):
    return (
        f"https://www.sdklab-biophysics-dzu.net/EVLncRNAs3/#/detail?id={record['ID']}"
    )


def description(record):
    return f"{record['Species']} {record['Class']}"


def species(record):
    return record["Species"]


def common_name(record):
    return record["Species"]


def disease_annotations(record, disease_mapping):
    """
    `disease_mapping` is normalised text -> (term id, match type), or None to
    attach nothing. A disease the mapping does not cover is an error: the
    mapping was built from another download, and skipping would lose data.
    """
    if disease_mapping is None:
        return []
    annotations = []
    for text, pmids in record["diseases"]:
        key = search.normalise(text)
        if key not in disease_mapping:
            raise ValueError(
                f"Disease {text!r} for {record['ID']} is not in the disease mapping"
            )
        term_id, match_type = disease_mapping[key]
        annotations.append(
            DiseaseAnnotation(
                disease_text=text,
                ontology_term_id=term_id,
                match_type=match_type,
                pubmed_ids=pmids,
                # Entry.database_name, i.e. rnc_database.descr
                assigned_by="EVLNCRNAS",
            )
        )
    return annotations


def as_entry(record, disease_mapping=None):
    return Entry(
        primary_id=primary_id(record),
        accession=accession(record),
        ncbi_tax_id=taxid(record),
        database="EVlncRNAs",
        sequence=sequence(record),
        regions=region_builder(record),
        rna_type=rna_type(record),
        url=url(record),
        seq_version="1",
        description=description(record),
        references=record["publications"],
        gene_synonyms=aliases(record),
        disease_annotations=disease_annotations(record, disease_mapping),
    )
