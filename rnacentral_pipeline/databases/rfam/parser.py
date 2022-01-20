# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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
import typing as ty
from pathlib import Path

from Bio import SeqIO

from ..data import Entry
from . import helpers, utils

# {
#     "rfam_acc": "RF00162",
#     "rfam_id": "SAM",
#     "type": "Cis-reg; riboswitch;",
#     "sequence_type": "full",
#     "description": "SAM riboswitch (S box leader)",
#     "ncbi_id": "69",
#     "rfamseq_acc": "CP013140.1",
#     "seq_start": "1228723",
#     "seq_end": "1228857",
#     "dbxrefs": "SO:0000035,GO:0010468",
#     "PMIDS": "PMID:10094622,PMID:12910260,PMID:12702767,PMID:15215334,PMID:16810258",
#     "version": "1",
#     "species": "Lysobacter enzymogenes",
#     "tax_string": "Bacteria; Proteobacteria; Gammaproteobacteria; Xanthomonadales; Xanthomonadaceae; Lysobacter.",
# }

# {
#     "id": "RF00162",
#     "name": "SAM",
#     "pretty_name": "SAM riboswitch (S box leader)",
#     "so_terms": "SO:0000035",
#     "rna_type": "Cis-reg; riboswitch;",
#     "description": "The SAM riboswitch is found upstream of a number of genes which code for proteins involved in methionine or cysteine biosynthesis in Gram-positive bacteria. The SAM riboswitch acts at the level of transcription termination control. The predicted structure consists of a complex stem-loop region followed by a single stem loop terminator region. An alternative and mutually exclusive form involves bases in the 3' segment of helix 1 with those in the 5' region of helix 5 to form a structure termed the anti-terminator form [1]. This entry represents the conserved core of the complex structure involving helices 1-4.",
#     "seed_count": "456",
#     "full_count": "6596",
#     "clan_id": "CL00012",
#     "length": "108",
# }


def as_entry(
    families: ty.Dict[str, ty.Dict[str, str]], sequences, data: ty.Dict[str, str]
) -> Entry:
    """
    Turn an entry from the JSON file into a Entry object for writing.
    """

    family = families[data["rfam_acc"]]
    return Entry(
        primary_id=helpers.primary_id(data),
        accession=helpers.accession(data),
        ncbi_tax_id=helpers.taxid(data),
        database="RFAM",
        sequence=helpers.sequence(sequences, data),
        regions=[],
        rna_type=helpers.rna_type(family),
        url=helpers.url(family),
        note_data=helpers.note(data),
        species=helpers.species(data),
        lineage=helpers.lineage(data),
        optional_id=helpers.optional_id(family),
        product=helpers.product(family),
        parent_accession=helpers.parent_accession(data),
        project="RFAM",
        experiment=helpers.experiment(data),
        description=helpers.description(family, data),
        mol_type=helpers.mol_type(data),
        seq_version=helpers.seq_version(data),
        is_composite="N",
        location_start=helpers.location_start(data),
        location_end=helpers.location_end(data),
        references=helpers.references(data),
    )


def load_mapping(handle: ty.TextIO) -> ty.Dict[str, ty.Dict[str, str]]:
    reader = csv.DictReader(handle, delimiter="\t")
    data = {}
    for row in reader:
        data[row["id"]] = row
    return data


def parse(
    family_file: ty.TextIO, sequence_info: ty.TextIO, fasta: Path
) -> ty.Iterable[Entry]:
    """
    Parse the JSON file of Rfam data and produce a generator of all Entry
    objects in the file.
    """

    indexed = SeqIO.index(str(fasta), "fasta")
    families = load_mapping(family_file)
    reader = csv.DictReader(sequence_info, delimiter="\t")
    for row in reader:
        yield as_entry(families, indexed, row)
