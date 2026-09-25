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

import attr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.helpers.publications import reference
from rnacentral_pipeline.databases.rfam import utils
from rnacentral_pipeline.databases.rfam.parser import as_entry, parse


def test_it_labels_y_rna_correctly():
    families = {
        "RF00019": {
            "id": "RF00019",
            "name": "Y_RNA",
            "pretty_name": "Y RNA",
            "so_terms": "SO:0000405",
        }
    }
    data = {
        "rfam_acc": "RF00019",
        "sequence_type": "seed",
        "description": "Y RNA",
        "ncbi_id": "7868",
        "rfamseq_acc": "AAVX01633839",
        "seq_start": "220",
        "seq_end": "328",
        "dbxrefs": "SO:0000405",
        "PMIDS": "10606662,7489501,10766734,18087752",
        "version": "1",
        "species": "Callorhinchus milii",
        "tax_string": (
            "Eukaryota Metazoa Chordata Craniata Vertebrata Chondrichthyes"
            " Holocephali Chimaeriformes Callorhinchidae Callorhinchus."
        ),
    }
    sequences = {
        "AAVX01633839/220-328": SeqRecord(
            Seq(
                "GGCCGGTCCGATGGTAGTGGGTTATCGTTGATATTTGCTTACAGAGTCAGTTACAG"
                "ATTTCCTTGTTCTCTCTTCCCCCCTTCTCACTGCTTCACTTGACTGGTCCTTT"
            )
        )
    }

    assert as_entry(families, sequences, data).rna_type == "SO:0000405"


def test_it_parses_all_data():
    with open("data/rfam/families.tsv") as mapping_handle, open(
        "data/rfam/sequence_info.tsv"
    ) as handle:
        assert (
            len(list(parse(mapping_handle, handle, "data/rfam/sequences.fasta")))
            == 1000
        )


def test_it_builds_first_entry_correctly():
    with open("data/rfam/families.tsv") as mapping_handle, open(
        "data/rfam/sequence_info.tsv"
    ) as handle:
        assert attr.asdict(
            next(parse(mapping_handle, handle, "data/rfam/sequences.fasta"))
        ) == attr.asdict(
            Entry(
                primary_id="RF00005",
                accession="KK113858.1:230594..230666:rfam",
                ncbi_tax_id=407821,
                database="RFAM",
                sequence=(
                    "GCCTTCGTGGTGTAGTGGTCAGCACACTTGACGCGTAACCGAGAGGTCCGTGGTTCGATTCTCG"
                    "GTGAAGGTG"
                ),
                regions=[],
                rna_type="SO:0000253",
                url="http://rfam.org/family/RF00005",
                note_data={
                    "Alignment": "full",
                    "SO": ["SO:0000253"],
                    "GO": ["GO:0030533"],
                },
                species="Stegodyphus mimosarum",
                lineage=(
                    "Eukaryota; Metazoa; Arthropoda; Chelicerata; Arachnida; Araneae; "
                    "Araneomorphae; Entelegynae; Eresoidea; Eresidae; Stegodyphus; "
                    "Stegodyphus mimosarum"
                ),
                references=[
                    reference(29112718),
                    reference(8256282),
                    reference(9023104),
                ],
                common_name=None,
                optional_id="tRNA",
                parent_accession="KK113858",
                description="Stegodyphus mimosarum tRNA",
                mol_type="full",
                location_start=230594,
                location_end=230666,
                seq_version="1",
                product="tRNA",
            )
        )


def test_it_can_parse_with_similar_accessions_correctly():
    with open("data/rfam/families.tsv") as mapping_handle, open(
        "data/rfam/sequence_info-duplicates.tsv"
    ) as handle:
        data = list(
            parse(mapping_handle, handle, "data/rfam/sequences-duplicates.fasta")
        )
    assert len(data) == 2
    assert data[0].accession == "CM000677.2:93286238..93286321:rfam"
    assert data[1].accession == "CM000677.2:93286321..93286238:rfam"
