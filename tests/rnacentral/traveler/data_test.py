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

from pathlib import Path

from rnacentral_pipeline.rnacentral.traveler import data


def test_paths_creates_expected_crw_paths():
    base = Path('data/traveler/crw/')
    path = data.TravelerPaths('URS00000F9D45_9606', 'd.5.e.H.sapiens.2', data.Source.crw, base)
    assert path.svg == Path('data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2.colored.svg')
    assert path.fasta == Path('data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2.fasta')
    assert path.overlaps == Path('data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2.overlaps')
    assert path.svg.exists()
    assert path.fasta.exists()
    assert path.overlaps.exists()
    # assert path.stk == Path('data/traveler/crw/URS00000F9D45_9606-d.5.e.H.sapiens.2.stk')


def test_paths_creates_expected_ribovision_paths():
    base = Path('data/traveler/ribovision/')
    path = data.TravelerPaths('URS0000AF4DA0', 'EC_LSU_3D', data.Source.ribovision, base)
    assert path.svg == Path('data/traveler/ribovision/URS0000AF4DA0-EC_LSU_3D.colored.svg')
    assert path.fasta == Path('data/traveler/ribovision/URS0000AF4DA0-EC_LSU_3D.fasta')
    assert path.overlaps == Path('data/traveler/ribovision/URS0000AF4DA0-EC_LSU_3D.overlaps')
    assert path.svg.exists()
    assert path.fasta.exists()
    assert path.overlaps.exists()
    # assert path.stk == Path('data/traveler/ribovision/URS0000AF4DA0-EC_LSU_3D.stk')


def test_paths_creates_expected_rfam_paths():
    base = Path('data/traveler/rfam/RF00162/')
    path = data.TravelerPaths('URS0000A7635A', 'RF00162', data.Source.rfam, base)
    assert path.svg == Path('data/traveler/rfam/RF00162/URS0000A7635A.colored.svg')
    assert path.fasta == Path('data/traveler/rfam/RF00162/URS0000A7635A.fasta')
    assert path.overlaps == Path('data/traveler/rfam/RF00162/URS0000A7635A.overlaps')
    assert path.svg.exists()
    assert path.fasta.exists()
    assert path.overlaps.exists()
    # assert path.stk == Path('data/traveler/rfam/RF00162/URS0000AF4DA0.stk')


def test_paths_creates_expected_gtrnadb_paths():
    base = Path('data/traveler/gtrnadb/')
    path = data.TravelerPaths('URS0000A0BF23', 'E-Gln', data.Source.gtrnadb, base)
    assert path.svg == Path('data/traveler/gtrnadb/URS0000A0BF23-E-Gln.colored.svg')
    assert path.fasta == Path('data/traveler/gtrnadb/URS0000A0BF23.fasta')
    assert path.overlaps == Path('data/traveler/gtrnadb/URS0000A0BF23-E-Gln.overlaps')
    assert path.svg.exists()
    assert path.fasta.exists()
    assert path.overlaps.exists()
    # assert path.stk == Path('data/traveler/gtrnadb/URS0000A0BF23-E-Gln.stk')
