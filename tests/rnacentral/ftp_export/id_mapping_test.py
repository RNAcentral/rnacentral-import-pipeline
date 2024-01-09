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

import pytest

from rnacentral_pipeline.rnacentral.ftp_export import id_mapping as ids
from tests.helpers import run_with_upi_taxid_constraint


@pytest.mark.parametrize(
    "data,expected",
    [
        (
            {
                "database": "ENA",
                "rna_type": "misc_RNA",
                "external_id": None,
                "optional_id": None,
                "accession": "EU334494.1:47..232:misc_RNA",
            },
            "EU334494.1:47..232:misc_RNA",
        ),
        (
            {
                "database": "ENSEMBL",
                "rna_type": "miRNA",
                "external_id": "ENSSSCT00000052986",
                "optional_id": "ENSSSCG00000037543.1",
                "accession": "ENSSSCT00000052986.1",
            },
            "ENSSSCT00000052986",
        ),
        (
            {
                "database": "HGNC",
                "rna_type": "lncRNA",
                "external_id": "A",
                "optional_id": "1",
                "accession": "HGNC:A",
            },
            "HGNC:A",
        ),
        (
            {
                "database": "PDBE",
                "rna_type": "rRNA",
                "external_id": "1S72",
                "optional_id": "1",
                "accession": "1S72_1_A",
            },
            "1S72_1",
        ),
    ],
)
def test_can_create_accession(data, expected):
    assert ids.accession(data) == expected


@pytest.mark.parametrize(
    "data,expected",
    [
        (
            {
                "gene": "A",
                "database": "PDBE",
                "rna_type": "rRNA",
            },
            "A",
        ),
        (
            {
                "gene": "something or other",
                "database": "PDBE",
                "rna_type": "rRNA",
            },
            "something or other",
        ),
        (
            {
                "gene": "first\tsecond",
                "database": "PDBE",
                "rna_type": "rRNA",
            },
            "first second",
        ),
        (
            {
                "gene": "piRNA-gene",
                "database": "ENA",
                "rna_type": "piRNA",
                "product": "piR-1",
            },
            "piR-1",
        ),
        (
            {
                "gene": "piR-1",
                "database": "ENSEMBL",
                "optional_id": "other",
                "rna_type": "piRNA",
            },
            "other",
        ),
        (
            {
                "gene": "gene1",
                "database": "ENSEMBL",
                "optional_id": "bob",
                "rna_type": "rRNA",
            },
            "bob",
        ),
        (
            {
                "gene": "gene1",
                "database": "MIRBASE",
                "optional_id": "hsa-mir-1",
                "rna_type": "miRNA",
            },
            "hsa-mir-1",
        ),
    ],
)
def test_can_generate_gene(data, expected):
    assert ids.gene(data) == expected


@pytest.mark.parametrize("name,expected", [("PDBE", "PDB"), ("HGNC", "HGNC")])
def test_can_create_databases(name, expected):
    assert ids.database({"database": name}) == expected


def test_as_entry_works_correctly():
    raw = {
        "upi": "a",
        "database": "PDBE",
        "external_id": "1S72",
        "optional_id": "1",
        "taxid": 102,
        "rna_type": "rRNA",
        "gene": None,
    }
    assert ids.as_entry(raw) == [
        "a",
        "PDB",
        "1S72_1",
        102,
        "rRNA",
        "",
    ]


@pytest.mark.parametrize(
    "rna_id,expected",
    [
        (
            "URS000018F875_9606",
            [
                [
                    "URS000018F875",
                    "ENSEMBL",
                    "ENST00000523510",
                    9606,
                    "lncRNA",
                    "ENSG00000254166.3",
                ],
                [
                    "URS000018F875",
                    "GENCODE",
                    "ENST00000523510",
                    9606,
                    "lncRNA",
                    "CASC19",
                ],
                ["URS000018F875", "HGNC", "HGNC:45089", 9606, "lncRNA", "PCAT2"],
                [
                    "URS000018F875",
                    "LNCIPEDIA",
                    "lnc-FAM84B-15:133",
                    9606,
                    "lncRNA",
                    "lnc-FAM84B-15",
                ],
                [
                    "URS000018F875",
                    "NONCODE",
                    "NONHSAT129016.2",
                    9606,
                    "lncRNA",
                    "NONHSAG051247.2",
                ],
                ["URS000018F875", "REFSEQ", "NR_119373", 9606, "lncRNA", "PCAT2"],
                # ['URS000018F875', 'GENCODE', 'ENST00000523510', 9606, 'lncRNA', 'ENSG00000254166.2'],
            ],
        ),
        (
            "URS0000672F0E_7955",
            [
                [
                    "URS0000672F0E",
                    "ENSEMBL",
                    "ENSDART00000171022",
                    7955,
                    "Y_RNA",
                    "ENSDARG00000100903.1",
                ],
                ["URS0000672F0E", "RFAM", "RF00019", 7955, "Y_RNA", ""],
            ],
        ),
        (
            "URS0000000096_9606",
            [
                [
                    "URS0000000096",
                    "ENA",
                    "DQ583192.1:1..30:ncRNA",
                    9606,
                    "piRNA",
                    "piR-50304",
                ],
            ],
        ),
        # (
        #     "URS000069C337_9606",
        #     [
        #         ["URS000069C337_9606", "ENSEMBL", "ENST00000401212", 9606, "pre-miRNA", "MIR298"],
        #         ["URS000069C337_9606", "GENECARDS", "", 9606, "pre-miRNA", "MIR298"],
        #         ["URS000069C337_9606", "MALACARDS", "", 9606, "pre-miRNA", "MIR298"],
        #         ["URS000069C337_9606", "MIRBASE", "", 9606, "pre-miRNA", "MIR298"],
        #         ["URS000069C337_9606", "REFSEQ", "", 9606, "pre-miRNA", "MIR298"],
        #     ]
        # ),
    ],
)
def test_can_create_expected_exports(rna_id, expected):
    entries = run_with_upi_taxid_constraint(
        rna_id, "files/ftp-export/id-mapping/id_mapping.sql", take_all=True
    )
    entries = sorted(ids.as_entry(e) for e in entries)
    assert entries == expected
