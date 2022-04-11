# -*- coding: utf-8 -*-

# pylint: disable=no-member

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

import os
import json
from collections import Counter

import attr
import pytest

from rnacentral_pipeline.databases.ensembl.metadata import assemblies as assem


@pytest.fixture(scope="module")
def assemblies():
    with open("config/databases.json", "r") as conn, open(
        "files/import-data/ensembl/assemblies.sql", "r"
    ) as query, open(
        "files/import-data/ensembl/example-locations.json", "r"
    ) as exs, open(
        "files/import-data/ensembl/known-assemblies.sql", "r"
    ) as kn:
        examples = json.load(exs)
        known = assem.load_known(os.environ["PGDATABASE"], kn)
        return list(assem.fetch(conn, query, examples, known))


def assembly_for(assemblies, taxid):
    found = [a for a in assemblies if a.taxid == taxid]
    assert len(found) == 1
    return found[0]


@pytest.mark.db
def test_it_can_load_known():
    with open("files/import-data/ensembl/known-assemblies.sql", "r") as kn:
        known = assem.load_known(os.environ["PGDATABASE"], kn)
    assert len(known) >= 440
    assert len(known[9606]) == 1
    assert attr.asdict(known[9606][0]) == attr.asdict(
        assem.AssemblyInfo(
            assembly_id="GRCh38",
            assembly_full_name="GRCh38.p13",
            gca_accession="GCA_000001405.28",
            assembly_ucsc="hg38",
            common_name="human",
            taxid=9606,
            ensembl_url="homo_sapiens",
            division="EnsemblVertebrates",
            blat_mapping=True,
            example=assem.AssemblyExample(chromosome="X", start=73819307, end=73856333),
        )
    )


@pytest.mark.db
def test_it_builds_a_valid_assembly(assemblies):
    val = assembly_for(assemblies, 9606)
    assert attr.asdict(val) == attr.asdict(
        assem.AssemblyInfo(
            assembly_id="GRCh38",
            assembly_full_name="GRCh38.p13",
            gca_accession="GCA_000001405.28",
            assembly_ucsc="hg38",
            common_name="human",
            taxid=9606,
            ensembl_url="homo_sapiens",
            division="EnsemblVertebrates",
            blat_mapping=True,
            example=assem.AssemblyExample(chromosome="X", start=73819307, end=73856333),
        )
    )
    assert val.subdomain == "ensembl.org"


@pytest.mark.db
def test_it_has_no_bacterial_assemblies(assemblies):
    divisions = set(a.division for a in assemblies)
    assert "EnsemblBacteria" not in divisions


@pytest.mark.db
@pytest.mark.parametrize(
    "taxid,count",
    [
        (4932, 0),
        (5127, 0),
        (546991, 1),
        (559292, 1),
        (6239, 1),
        (6669, 1),
        (7227, 1),
        (8090, 1),
    ],
)
def test_it_has_one_assembly_per_taxid(assemblies, taxid, count):
    val = [a for a in assemblies if a.taxid == taxid]
    assert len(val) == count, "Issue with counts for %i" % taxid


@pytest.mark.db
@pytest.mark.parametrize("key", ["taxid", "assembly_id"])
def test_it_never_has_more_than_one_assembly_unique_item(assemblies, key):
    counts = Counter(getattr(a, key) for a in assemblies)
    max_item = counts.most_common(n=1)[0]
    assert max_item[1] == 1, "Too many counts for: %s, %i" % max_item


@pytest.mark.db
@pytest.mark.parametrize(
    "taxid,division,assembly_id",
    [
        (546991, "EnsemblFungi", "EF2"),
        (559292, "EnsemblFungi", "R64-1-1"),
        (6239, "EnsemblVertebrates", "WBcel235"),
        (6669, "EnsemblMetazoa", "V1.0"),
        (7227, "EnsemblVertebrates", "BDGP6.28"),
        (8090, "EnsemblVertebrates", "ASM223467v1"),
    ],
)
def test_it_uses_correct_sources_for_duplicates(
    assemblies, taxid, division, assembly_id
):
    val = assembly_for(assemblies, taxid)
    assert val.division == division
    assert val.assembly_id == assembly_id
