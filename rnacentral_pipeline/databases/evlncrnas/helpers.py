import numpy as np

from rnacentral_pipeline.databases.data import Entry, Exon, SequenceRegion
from rnacentral_pipeline.databases.helpers import phylogeny as phy


def primary_id(record):
    return "EVLNCRNAS:" + record["ID"]


def accession(record):
    return "EVLNCRNAS:" + record["ID"]


def taxid(record):
    return int(record["taxid"])


def sequence(record):
    return record["sequence"].upper().replace("U", "T")


def aliases(record):
    print(f'Name: {record["external_id"]}, aliases: {record["Aliases"]}')

    if record["Aliases"] is None:
        return [str(record["external_id"])]

    aliases = [str(record["external_id"])]
    aliases.extend(str(record["Aliases"]).split(","))

    return aliases


def region_builder(info):
    if not info["region_start"] or not info["region_stop"]:
        return []
    return [
        SequenceRegion(
            chromosome=info["chromosome"],
            strand=info["Chain"],
            exons=[
                Exon(start=int(info["region_start"]), stop=int(info["region_stop"]))
            ],
            assembly_id=info["assembly_id"],
            coordinate_system="1-start, fully-closed",
        )
    ]


def rna_type(record):
    return record["Class"]
    pass


def url(record):
    return (
        "https://www.sdklab-biophysics-dzu.net/EVLncRNAs2/index.php/Home/Browsc/rna.html?id="
        + record["ID"]
    )


def description(record):
    return f"{record['Species']} {record['Class']} from EVlncRNAs"


def species(record):
    return record["Species"]


def common_name(record):
    return record["Species"]


def lineage(record):
    return phy.lineage(int(record["taxid"]))


def as_entry(record):
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
        species=species(record),
        common_name=common_name(record),
        lineage=lineage(record),
        references=record["publications"],
        gene_synonyms=aliases(record),
    )
