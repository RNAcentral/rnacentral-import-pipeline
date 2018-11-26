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

from collections import Counter

import attr
import pytest

from rnacentral_pipeline.databases.ensembl.metadata import assemblies as assem


@pytest.fixture(scope='module')
def assemblies():
    with open('config/databases.json', 'r') as conn, \
            open('files/import-data/ensembl/assemblies.sql', 'r') as query:
        return list(assem.fetch(conn, query, assem.EXAMPLE_LOCATIONS))


def assembly_for(assemblies, taxid):
    found = [a for a in assemblies if a.taxid == taxid]
    assert len(found) == 1
    return found[0]


def test_it_builds_a_valid_assembly(assemblies):
    val = assembly_for(assemblies, 9606)
    assert attr.asdict(val) == attr.asdict(assem.AssemblyInfo(
        assembly_id='GRCh38',
        assembly_full_name='GRCh38.p12',
        gca_accession='GCA_000001405.27',
        assembly_ucsc='hg38',
        common_name='human',
        taxid=9606,
        ensembl_url='homo_sapiens',
        division='Ensembl',
        blat_mapping=True,
        example=assem.AssemblyExample(
            chromosome='X',
            start=73819307,
            end=73856333,
        ),
    ))
    assert val.subdomain == 'ensembl.org'


def test_it_has_no_bacterial_assemblies(assemblies):
    divisions = set(a.division for a in assemblies)
    assert 'EnsemblBacteria' not in divisions


def test_it_has_one_assembly_per_taxid(assemblies):
    counts = Counter(a.taxid for a in assemblies)
    max_item = counts.most_common(n=1)[0]
    assert max_item[1] == 1
    for taxid in [7227, 6239, 8090]:
        assert counts[taxid] == 1


def test_it_has_one_assembly_per_assembly_id(assemblies):
    counts = Counter(a.assembly_id for a in assemblies)
    max_item = counts.most_common(n=1)[0]
    assert max_item[1] == 1


@pytest.mark.parametrize('taxid', [  # pylint: disable=no-member
    (4932),
    (5127),
])
def test_it_does_not_have_unexpected_taxids(assemblies, taxid):
    assert len([a for a in assemblies if a.taxid == taxid]) == 0


@pytest.mark.parametrize('taxid,division,assembly_id', [
    (7227, 'Ensembl', 'BDGP6'),
    (6239, 'Ensembl', 'WBcel235'),
    (559292, 'EnsemblFungi', 'R64-1-1'),
    (6669, 'EnsemblMetazoa', 'V1.0'),
    (546991, 'EnsemblFungi', 'EF2'),
    (8090, 'Ensembl', 'ASM223467v1'),
])
def test_it_uses_correct_sources_for_duplicates(assemblies, taxid, division,
                                                assembly_id):
    val = assembly_for(assemblies, taxid)
    assert val.division == division
    assert val.assembly_id == assembly_id
