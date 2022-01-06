# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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
from pathlib import Path
import datetime as dt

import pytest

from rnacentral_pipeline.databases.psicquic import parser
from rnacentral_pipeline.databases.data import (
    Entry,
    Interaction,
    InteractionIdentifier,
    Interactor,
)
from rnacentral_pipeline.databases.helpers import publications as pub


@pytest.fixture(scope="module")
def data():
    path = Path("data/psicquic/data.tsv")
    return list(parser.parse(path, os.environ["PGDATABASE"]))


def test_can_parse_all_entries(data):
    assert len(data) == 9


def test_can_parse_correctly(data):
    assert data[0] == Entry(
        primary_id="PSICQUIC:URS000000B1C9_9606",
        accession="PSICQUIC:URS000000B1C9_9606",
        ncbi_tax_id=9606,
        database="PSICQUIC",
        sequence="TGAGGTAGGAGGTTGTATAGTT",
        regions=[],
        rna_type="miRNA",
        url="http://www.ebi.ac.uk/Tools/webservices/psicquic/view/main.xhtml",
        seq_version="1",
        description="Homo sapiens (human) hsa-let-7e-5p",
        references=[
            pub.reference("PMID:23671334"),
            pub.reference("PMID:30670152"),
        ],
        species="Homo sapiens",
        common_name="human",
        lineage="Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi; Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini; Catarrhini; Hominidae; Homo; Homo sapiens",
        interactions=[
            Interaction(
                ids=(
                    InteractionIdentifier(
                        key="psicquic", value="URS000000B1C9_9606-0", name=""
                    ),
                ),
                interactor1=Interactor(
                    id=InteractionIdentifier(
                        key="RNAcentral", value="URS000000B1C9_9606", name=None
                    ),
                    alt_ids=(),
                    aliases=(
                        InteractionIdentifier(
                            key="RNAcentral",
                            value="Homo sapiens (human) hsa-let-7e-5p",
                            name="recommended name",
                        ),
                    ),
                    taxid=9606,
                    biological_role=(
                        InteractionIdentifier(
                            key="psi-mi", value="MI:0499", name="unspecified role"
                        ),
                    ),
                    experimental_role=(
                        InteractionIdentifier(
                            key="psi-mi", value="MI:0499", name="unspecified role"
                        ),
                    ),
                    interactor_type=(
                        InteractionIdentifier(
                            key="psi-mi", value="MI:0320", name="ribonucleic acid"
                        ),
                    ),
                    xrefs=(
                        InteractionIdentifier(
                            key="go",
                            value="GO:0035278",
                            name="miRNA mediated inhibition of translation",
                        ),
                    ),
                    annotations="-",
                    features=(),
                    stoichiometry=None,
                    participant_identification=(),
                ),
                interactor2=Interactor(
                    id=InteractionIdentifier(
                        key="uniprotkb", value="P40763", name=None
                    ),
                    alt_ids=(),
                    aliases=(),
                    taxid=-3,
                    biological_role=(
                        InteractionIdentifier(
                            key="psi-mi", value="MI:0499", name="unspecified role"
                        ),
                    ),
                    experimental_role=(
                        InteractionIdentifier(
                            key="psi-mi", value="MI:0499", name="unspecified role"
                        ),
                    ),
                    interactor_type=(
                        InteractionIdentifier(
                            key="psi-mi", value="MI:0329", name="unknown participant"
                        ),
                    ),
                    xrefs=(),
                    annotations="-",
                    features=(),
                    stoichiometry=None,
                    participant_identification=(),
                ),
                methods=(),
                types=(
                    InteractionIdentifier(
                        key="psi-mi", value="MI:0915", name="physical association"
                    ),
                ),
                xrefs=(),
                annotations=(),
                confidence=(),
                source_database=(
                    InteractionIdentifier(
                        key="psi-mi", value="MI:2320", name="aruk-ucl"
                    ),
                ),
                is_negative=True,
                publications=(pub.reference("PMID:30670152"),),
                create_date=None,
                update_date=dt.date(2020, 1, 27),
                host_organisms=None,
            )
        ],
    )
