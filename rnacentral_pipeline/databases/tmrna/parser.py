# -*- coding: utf-8 -*-

"""
Copyright [2009-${2024}] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import csv
import sys
import typing as ty

from Bio import Entrez

from rnacentral_pipeline.databases.data import Entry, SequenceFeature
from rnacentral_pipeline.databases.helpers import gtdb

csv.field_size_limit(sys.maxsize)  # It contains huge rows
Entrez.email = "rnacentral@gmail.com"

name_mapping = {
    "Acceptor": "tmrna_acceptor",
    "Body": "tmrna_body",
    "CCAequiv": "tmrna_ccaequiv",
    "Coding": "tmrna_coding_region",
    "Exon1": "tmrna_exon",
    "Exon2": "tmrna_exon",
    "GpI": "tmrna_gpi",
    "IVS": "tmrna_ivs",
    "TagCDS": "tmrna_tagcds",
}

NO_CODING_SEQUENCE = ["frameshift", "undetermined"]

ORGANELLE = ["Chromatophores", "Mitochondria", "Plastids"]


def inferred_species(id: str) -> str:
    raw = id.split(".")[0]
    return raw.replace("__", " ")


def features(raw: ty.Dict[str, str]) -> ty.List[SequenceFeature]:
    features = []
    segements = raw["Segments"].split(",")
    for segement in segements:
        name, range = segement.split(":")
        start, stop = range.split("-")
        # Convert to zero based
        location = [int(start) - 1, int(stop)]
        metadata = {}
        if name == "TagCDS":
            if raw["Tag"] in NO_CODING_SEQUENCE:
                metadata["orf_summary"] = raw["Tag"]
            else:
                coding = raw["Tag"]
                metadata["orf_summary"] = "Coding sequence"
                metadata["has_stop"] = coding[-1] == "*"
                if coding[-1] == "*":
                    coding = coding[0:-1]
                metadata["coding_sequence"] = coding
        feature = SequenceFeature(
            name=name_mapping[name],
            feature_type=name_mapping[name],
            location=location,
            sequence="",
            provider="TMRNA_WEB",
            metadata=metadata,
        )
        features.append(feature)
    return features


def parse(raw: ty.IO) -> ty.Iterable[Entry]:
    reader = csv.DictReader(raw, delimiter="\t")
    rows = list(reader)
    for row in rows:
        accessions = row["Instances"].split(",")
        assert len(accessions) == int(row["InstanceCt"])
        species = inferred_species(row["#ID"])
        tax_string = ",".join([row["Taxonomy"], species])
        tax_id = gtdb.phylogeny_to_taxid(tax_string)
        for accession in accessions:
            accession_parts = accession.split("/")
            note = {
                "tmrna_form": row["Form"],
            }
            yield Entry(
                accession=f"tmrna:{accession}",
                primary_id=f"tmrna:{accession}",
                ncbi_tax_id=tax_id,
                database="TMRNA_WEB",
                sequence=row["Sequence"].upper(),
                regions=[],
                rna_type="SO:0000584",
                url="",
                seq_version="1",
                note_data=note,
                parent_accession=accession_parts[0],
                features=features(row),
                inference=row["Taxonomy"],
                description=f"{species} {row['Form']} tmRNA",
            )
