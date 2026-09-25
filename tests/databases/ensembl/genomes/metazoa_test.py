# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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
import pytest

from rnacentral_pipeline.databases import data as dat
from rnacentral_pipeline.databases.ensembl import metazoa
from rnacentral_pipeline.databases.helpers import publications as pubs

from . import helpers

# Amel_4.5 (this fixture's original assembly) was replaced by Amel_HAv3.1,
# with an entirely different gene id scheme (ENSRNA0227... -> RefSeq-style
# ids) - re-picked five genes covering different ncRNA types from the
# current assembly rather than trying to track the ~30 original ids.


@pytest.fixture(scope="module")  # pylint: disable=no-member
def apis_1():
    return helpers.parse(
        metazoa.parse,
        "test-data/ensembl_metazoa/Apis_mellifera.Amel_HAv3.1.primary_assembly.CM009931.2.dat",
    )


@pytest.mark.parametrize(
    "filename,count",
    [
        (
            "test-data/ensembl_metazoa/Apis_mellifera.Amel_HAv3.1.primary_assembly.CM009931.2.dat",
            5,
        ),
    ],
)
def test_can_parse_all_entries(filename, count):
    assert len(helpers.parse(metazoa.parse, filename)) == count


def test_gets_all_ids(apis_1):
    assert set(d.primary_id for d in apis_1) == {
        "XR_003306550",  # snRNA
        "XR_003306551",  # snoRNA
        "XR_411947",  # lncRNA
        "GeneID_732512_t1",  # miRNA
        "GeneID_107966081_t1",  # tRNA
    }


@pytest.mark.parametrize(
    "accession,rna_type",
    [
        ("ENSEMBL_METAZOA:XR_003306550", "SO:0000673"),
        ("ENSEMBL_METAZOA:GeneID_107966081_t1", "SO:0000673"),
    ],
)
def test_can_assign_expected_rna_types(apis_1, accession, rna_type):
    assert helpers.entry_for(apis_1, accession).rna_type == rna_type


@pytest.mark.parametrize(
    "accession,description",
    [
        (
            "ENSEMBL_METAZOA:XR_003306550",
            "Apis mellifera (Honey bee) U1 spliceosomal RNA",
        ),
    ],
)
def test_can_get_expected_descriptions(apis_1, accession, description):
    assert helpers.entry_for(apis_1, accession).description == description


def test_can_get_expected_data(apis_1):
    assert attr.asdict(
        helpers.entry_for(apis_1, "ENSEMBL_METAZOA:XR_003306550")
    ) == attr.asdict(
        dat.Entry(
            primary_id="XR_003306550",
            accession="ENSEMBL_METAZOA:XR_003306550",
            ncbi_tax_id=7460,
            database="ENSEMBL_METAZOA",
            sequence=(
                "TAACTTACTTGGCGCGGAGGATACCGTGATCACGAAGGCGGTTCCTTCTGGGCGAGGCT"
                "CTTCCATTGCACTTAGGTAGAGCTGAACCTTGCGAATACTCCTAATGTGGGTATGTCGA"
                "GCGCACAATTTTTGGTAGTCGGGACCTGCGTTCGCGCTGTCCCGGA"
            ),
            regions=[
                dat.SequenceRegion(
                    chromosome="CM009931.2",
                    strand=1,
                    exons=[dat.Exon(start=31, stop=194)],
                    assembly_id="Amel_HAv3.1",
                    coordinate_system=dat.CoordinateSystem.one_based(),
                ),
            ],
            rna_type="SO:0000673",
            url="",
            seq_version="1",
            note_data={},
            xref_data={
                "GenBank_transcript": ["XR_003306550.1"],
                "RFAM_trans_name": ["RF00003"],
            },
            species="Apis mellifera",
            common_name="Honey bee",
            lineage=(
                "Eukaryota; Metazoa; Ecdysozoa; Arthropoda; Altocrustacea; "
                "Allotriocarida; Hexapoda; Insecta; Pterygota; Neoptera; "
                "Eumetabola; Endopterygota; Hymenoptera; Apocrita; Aculeata; "
                "Apoidea; Anthophila; Apidae; Apis; Apis mellifera"
            ),
            gene="LOC113219421",
            description="Apis mellifera (Honey bee) U1 spliceosomal RNA",
            references=[pubs.reference("doi:10.1093/nar/gkx1011")],
        )
    )


def test_does_not_create_ncRNA_rna_type(apis_1):
    for entry in apis_1:
        assert entry.rna_type != "ncRNA"
