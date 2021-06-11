# -*- coding: utf-8 -*-

"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

import os
from datetime import date
import collections as coll

import attr
import pytest

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.intact import parser
from rnacentral_pipeline.databases.helpers import publications as pubs


@pytest.fixture(scope="module")
def sample():
    with open("data/intact/sample.txt", "r") as raw:
        return list(parser.parse(raw, os.environ["PGDATABASE"]))


def test_can_parse_all_data(sample):
    assert len(sample) == 4


def test_creates_entries_with_expected_ids(sample):
    ids = sorted(e.primary_id for e in sample)
    assert ids == [
        "INTACT:URS0000077671_559292",
        "INTACT:URS0000182FAB_559292",
        "INTACT:URS000020D517_559292",
        "INTACT:URS00002AFD52_559292",
    ]


def test_correctly_groups_data(sample):
    val = {e.primary_id: len(e.interactions) for e in sample}
    assert val == {
        "INTACT:URS0000077671_559292": 1,
        "INTACT:URS0000182FAB_559292": 1,
        "INTACT:URS000020D517_559292": 1,
        "INTACT:URS00002AFD52_559292": 4,
    }


def test_produces_correct_data(sample):
    with open("data/intact/sample.txt", "r") as raw:
        val = next(parser.parse_interactions(raw))
    i1 = data.Interactor(
        id=data.InteractionIdentifier("intact", "EBI-10921362", None),
        alt_ids=[
            data.InteractionIdentifier("rnacentral", "URS00002AFD52_559292", None)
        ],
        aliases=[
            data.InteractionIdentifier("psi-mi", "snr18 yeast", "display_short"),
            data.InteractionIdentifier("psi-mi", "EBI-10921362", "display_long"),
        ],
        taxid=559292,
        biological_role=[
            data.InteractionIdentifier("psi-mi", "MI:0499", "unspecified role")
        ],
        experimental_role=[data.InteractionIdentifier("psi-mi", "MI:0498", "prey")],
        interactor_type=[
            data.InteractionIdentifier("psi-mi", "MI:0609", "small nucleolar rna")
        ],
        xrefs=[],
        annotations="-",
        features=[data.InteractionIdentifier("32p radiolabel", "?-?", None)],
        stoichiometry=None,
        participant_identification=[
            data.InteractionIdentifier("psi-mi", "MI:0396", "predetermined participant")
        ],
    )
    i2 = data.Interactor(
        id=data.InteractionIdentifier("uniprotkb", "P15646", None),
        alt_ids=[
            data.InteractionIdentifier("intact", "EBI-6838", None),
            data.InteractionIdentifier("uniprotkb", "P89890", None),
            data.InteractionIdentifier("uniprotkb", "D6VRX5", None),
        ],
        aliases=[
            data.InteractionIdentifier("psi-mi", "fbrl_yeast", "display_long"),
            data.InteractionIdentifier("uniprotkb", "NOP1", "gene name"),
            data.InteractionIdentifier("psi-mi", "NOP1", "display_short"),
            data.InteractionIdentifier("uniprotkb", "LOT3", "gene name synonym"),
            data.InteractionIdentifier("uniprotkb", "YDL014W", "locus name"),
            data.InteractionIdentifier("uniprotkb", "D2870", "orf name"),
            data.InteractionIdentifier(
                "uniprotkb",
                "U3 small nucleolar RNA-associated protein NOP1",
                "gene name synonym",
            ),
            data.InteractionIdentifier(
                "uniprotkb", "Histone-glutamine methyltransferase", "gene name synonym"
            ),
        ],
        taxid=559292,
        biological_role=[
            data.InteractionIdentifier("psi-mi", "MI:0499", "unspecified role")
        ],
        experimental_role=[data.InteractionIdentifier("psi-mi", "MI:0498", "prey")],
        interactor_type=[data.InteractionIdentifier("psi-mi", "MI:0326", "protein")],
        xrefs=[
            data.InteractionIdentifier(
                "go", "GO:0008649", "rRNA methyltransferase activity"
            ),
            data.InteractionIdentifier("go", "GO:0031428", "box C/D snoRNP complex"),
            data.InteractionIdentifier("go", "GO:0032040", "small-subunit processome"),
            data.InteractionIdentifier("go", "GO:0003723", "RNA binding"),
            data.InteractionIdentifier(
                "go", "GO:0000494", "box C/D snoRNA 3'-end processing"
            ),
            data.InteractionIdentifier("refseq", "NP_010270.1", None),
            data.InteractionIdentifier("sgd", "S000002172", None),
            data.InteractionIdentifier("interpro", "IPR000692", "Fibrillarin"),
            data.InteractionIdentifier("interpro", "IPR020813", None),
            data.InteractionIdentifier("rcsb pdb", "5WYJ", None),
            data.InteractionIdentifier("rcsb pdb", "5WYK", None),
            data.InteractionIdentifier(
                "go", "GO:0006356", "regulation of transcription by RNA polymerase I"
            ),
            data.InteractionIdentifier(
                "go", "GO:1990259", "histone-glutamine methyltransferase activity"
            ),
            data.InteractionIdentifier("go", "GO:0005730", "nucleolus"),
            data.InteractionIdentifier("go", "GO:0006364", "rRNA processing"),
            data.InteractionIdentifier("go", "GO:0031167", "rRNA methylation"),
            data.InteractionIdentifier("go", "GO:0043144", "snoRNA processing"),
            data.InteractionIdentifier(
                "go", "GO:1990258", "histone glutamine methylation"
            ),
            data.InteractionIdentifier("interpro", "IPR029063", None),
            data.InteractionIdentifier("mint", "P15646", None),
            data.InteractionIdentifier("go", "GO:0000451", "rRNA 2'-O-methylation"),
            data.InteractionIdentifier("go", "GO:0005654", "nucleoplasm"),
            data.InteractionIdentifier(
                "go", "GO:0008171", "O-methyltransferase activity"
            ),
            data.InteractionIdentifier("go", "GO:0015030", "Cajal body"),
            data.InteractionIdentifier("reactome", "R-SCE-6791226", None),
            data.InteractionIdentifier("dip", "DIP-698N", None),
            data.InteractionIdentifier("rcsb pdb", "5WLC", None),
            data.InteractionIdentifier("go", "GO:0030686", "90S preribosome"),
            data.InteractionIdentifier("rcsb pdb", "6ND4", None),
        ],
        annotations="crc64:56A8B958A7B6066E",
        features=[data.InteractionIdentifier("protein a tag", "n-n", None)],
        stoichiometry=None,
        participant_identification=[
            data.InteractionIdentifier("psi-mi", "MI:0396", "predetermined participant")
        ],
    )

    assert attr.asdict(val) == attr.asdict(
        data.Interaction(
            ids=[data.InteractionIdentifier("intact", "EBI-11665247", None)],
            interactor1=i1,
            interactor2=i2,
            methods=[],
            types=[
                data.InteractionIdentifier("psi-mi", "MI:0915", "physical association")
            ],
            xrefs=[],
            annotations=[
                data.InteractionIdentifier("figure legend", "Fig 2, Fig 3B", None)
            ],
            confidence=[data.InteractionIdentifier("intact-miscore", "0.74", None)],
            source_database=[data.InteractionIdentifier("psi-mi", "MI:0471", "MINT")],
            is_negative=False,
            publications=[pubs.reference(11726521)],
            create_date=date(2003, 7, 8),
            update_date=date(2016, 3, 23),
            host_organisms=-1,
        )
    )
