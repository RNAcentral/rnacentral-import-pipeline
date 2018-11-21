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


import attr

import pytest

from rnacentral_pipeline.rnacentral import genome_mapping as gm
from rnacentral_pipeline.databases.data import regions


def parse(assembly_id, filename):
    with open(filename, 'r') as raw:
        return list(gm.parse(assembly_id, raw))


def parse_of(assembly_id, filename, upi):
    return [h for h in parse(assembly_id, filename) if h.upi == upi]


def select_hits(assembly_id, filename):
    with open(filename, 'r') as raw:
        return list(gm.select_hits(assembly_id, raw))


def results_for(assembly_id, filename, upi):
    return [h for h in select_hits(assembly_id, filename) if h.upi == upi]


def result_for(assembly_id, filename, upi):
    results = results_for(assembly_id, filename, upi)
    assert len(results) == 1
    return results[0]


@pytest.mark.parametrize('assembly_id,filename,count', [
    ('human', 'data/genome-mapping/results.psl', 47195),
])
def test_parses_correct_data(assembly_id, filename, count):
    assert len(parse(assembly_id, filename)) == count


def test_gets_correct_upis():
    hits = select_hits('human', 'data/genome-mapping/results.psl')
    val = set(h.upi for h in hits)
    assert val == {
        'URS0000032237_9606',
        'URS0000034CE9_9606',
        'URS000007D899_9606',
        'URS00000B856F_9606',
        'URS00000CF258_9606',
        'URS000013CF05_9606',
        'URS0000169307_9606',
        'URS000019683B_9606',
        'URS0000233572_9606',
        'URS00002412DD_9606',
        'URS000024BA85_9606',
        'URS000024DA9C_9606',
        'URS000027D6D4_9606',
        'URS00002EF678_9606',
        'URS000031BFD7_9606',
        'URS000033979A_9606',
        'URS000034CE59_9606',
        'URS00003C95EB_9606',
        'URS00003D311C_9606',
        'URS00003ED72C_9606',
        'URS0000463768_9606',
        'URS00004AA950_9606',
        'URS000050824C_9606',
        'URS0000551482_9606',
        'URS0000567331_9606',
        'URS0000584A85_9606',
        'URS00005ADEFA_9606',
        'URS00005C4E95_9606',
        'URS00006D05AC_9606',
        'URS000073FCEB_9606',
        'URS000075EF8B_9606',
        'URS00007C3E04_9606',
        'URS00007CC811_9606',
        'URS00007CCDF8_9606',
        'URS00008C2236_9606',
        'URS00008C2810_9606',
        'URS0000907D4C_9606',
        'URS0000924419_9606',
        'URS000093C0D6_9606',
        'URS0000964436_9606',
        'URS000097430A_9606',
        'URS0000A7635E_9606',
        'URS0000B70A56_9606',
        'URS0000BB3005_9606',
        'URS0000BC164F_9606',
        'URS0000CB4B22_9606',
        'URS0000CD47EC_9606',
        'URS0000CE4269_9606',
        'URS0000CED72F_9606',
        'URS0000CF712B_9606',
        'URS0000CF843A_9606',
        'URS0000D051B0_9606',
        'URS0000D0D132_9606',
        'URS0000D35D47_9606',
        'URS0000D37792_9606',
        'URS0000D56C31_9606',
        'URS0000D56C43_9606',
        'URS0000D56CE0_9606',
        'URS0000D56E36_9606',
        'URS0000D56E4F_9606',
        'URS0000D57299_9606',
        'URS0000D572C8_9606',
        'URS0000D57454_9606',
        'URS0000D574F6_9606',
        'URS0000D57550_9606',
        'URS0000D57551_9606',
        'URS0000D57619_9606',
        'URS0000D57675_9606',
        'URS0000D57757_9606',
        'URS0000D577F3_9606',
        'URS0000D579AF_9606',
        'URS0000D57A7D_9606',
        'URS0000D57ADC_9606',
        'URS0000D57AF2_9606',
        'URS0000D57B22_9606',
        'URS0000D57C6C_9606',
        'URS0000D57E65_9606',
        'URS0000D57EBD_9606',
        'URS0000D57ECA_9606',
        'URS0000D582A5_9606',
        'URS0000D58337_9606',
        'URS0000D58361_9606',
        'URS0000D58375_9606',
        'URS0000D584E4_9606',
        'URS0000D58593_9606',
        'URS0000D585B5_9606',
        'URS0000D58715_9606',
        'URS0000D58726_9606',
        'URS0000D5889B_9606',
        'URS0000D588E0_9606',
        'URS0000D5896F_9606',
        'URS0000D58A24_9606',
        'URS0000D58B1C_9606',
        'URS0000D58CD2_9606',
        'URS0000D58D1E_9606',
        'URS0000D58D68_9606',
        'URS0000D58D98_9606',
        'URS0000D58FD5_9606',
        'URS0000D5900A_9606',
        'URS0000D59281_9606',
        'URS0000D59343_9606',
        'URS0000D593C8_9606',
        'URS0000D595A4_9606',
        'URS0000D59668_9606',
        'URS0000D59825_9606',
        'URS0000D5999A_9606',
        'URS0000D59D35_9606',
        'URS0000D59D36_9606',
        'URS0000D59D65_9606',
        'URS0000D59EA3_9606',
        'URS0000D59F7A_9606',
        'URS0000D5A24A_9606',
        'URS0000D5A284_9606',
        'URS0000D5A2E1_9606',
        'URS0000D5A3C6_9606',
        'URS0000D5A43C_9606',
        'URS0000D5A583_9606',
        'URS0000D5A6A4_9606',
        'URS0000D5A846_9606',
        'URS0000D5ACF7_9606',
        'URS0000D5AD75_9606',
        'URS0000D5AD9C_9606',
        'URS0000D5AEF3_9606',
        'URS0000D5AFF1_9606',
        'URS0000D5B144_9606',
        'URS0000D5B3D7_9606',
        'URS0000D5B428_9606',
        'URS0000D5B486_9606',
        'URS0000D5B4C1_9606',
        'URS0000D5B4F1_9606',
        'URS0000D5B517_9606',
        'URS0000D5B965_9606',
        'URS0000D5B989_9606',
        'URS0000D5BA30_9606',
        'URS0000D5BBA8_9606',
        'URS0000D5BBAA_9606',
        'URS0000D5BC6E_9606',
        'URS0000D5BDDB_9606',
        'URS0000D5BE08_9606',
        'URS0000D5BFF2_9606',
        'URS0000D5C169_9606',
        'URS0000D5C30A_9606',
        'URS0000D5C3F8_9606',
        'URS0000D5C43E_9606',
        'URS0000D5C543_9606',
        'URS0000D5C63A_9606',
        'URS0000D5C81B_9606',
        'URS0000D5C835_9606',
        'URS0000D5C84D_9606',
        'URS0000D5C8A0_9606',
        'URS0000D5C8E8_9606',
        'URS0000D5C94A_9606',
        'URS0000D5C966_9606',
        'URS0000D5C9D5_9606',
        'URS0000D5C9E4_9606',
        'URS0000D5CBA1_9606',
        'URS0000D5CDF6_9606',
        'URS0000D5D205_9606',
        'URS0000D5D2B3_9606',
        'URS0000D5D33F_9606',
        'URS0000D5D42C_9606',
        'URS0000D5D561_9606',
        'URS0000D5D669_9606',
        'URS0000D5D6C8_9606',
        'URS0000D5D6F7_9606',
        'URS0000D5D726_9606',
        'URS0000D5D797_9606',
        'URS0000D5D7A2_9606',
        'URS0000D5D7E9_9606',
        'URS0000D5D88F_9606',
        'URS0000D5D8C3_9606',
        'URS0000D5DB34_9606',
        'URS0000D5DBAE_9606',
        'URS0000D5DD92_9606',
        'URS0000D5DF9E_9606',
        'URS0000D5E072_9606',
        'URS0000D5E0AE_9606',
        'URS0000D5E20B_9606',
        'URS0000D5E32C_9606',
        'URS0000D5E35B_9606',
        'URS0000D77905_9606',
        'URS0000D77926_9606',
        'URS0000D7793C_9606',
        'URS0000D77F51_9606',
        'URS0000D77F56_9606',
        'URS0000D77FAA_9606',
        'URS0000D77FC6_9606',
        'URS0000D78002_9606',
        'URS0000D7802D_9606',
        'URS0000D78084_9606',
        'URS0000D780C8_9606',
        'URS0000D78103_9606',
        'URS0000D78123_9606',
        'URS0000D926BB_9606',
        'URS0000D9ADDE_9606',
        'URS0000D9D3FE_9606',
        'URS0000DA77B0_9606',
        'URS0000DB164A_9606',
        'URS0000DB2205_9606',
    }


def test_can_build_correct_data_structures():
    val = parse_of('human', 'data/genome-mapping/results.psl', 'URS0000032237_9606')
    assert len(val) == 1
    assert attr.asdict(val[0]) == attr.asdict(gm.Hit(
        assembly_id='human',
        chromosome='MT',
        upi='URS0000032237_9606',
        sequence_length=19,
        matches=19,
        target_insertions=0,
        strand=1,
        exons=[regions.Exon(start=15935, stop=15953)]
    ))

def test_can_build_correct_for_minus_strand():
    val = parse_of('human', 'data/genome-mapping/results.psl', 'URS00000CF258_9606')
    assert len(val) == 1
    assert attr.asdict(val[0]) == attr.asdict(gm.Hit(
        assembly_id='human',
        chromosome='Y',
        upi='URS00000CF258_9606',
        sequence_length=16,
        matches=16,
        target_insertions=0,
        strand=-1,
        exons=[regions.Exon(start=17624526, stop=17624541)]
    ))



def test_can_correctly_parse_with_several_exons():
    val = parse_of('human', 'data/genome-mapping/results.psl', 'URS000007D899_9606')
    assert len(val) == 29
    assert attr.asdict(val[0]) == attr.asdict(gm.Hit(
        assembly_id='human',
        chromosome='21',
        upi='URS000007D899_9606',
        sequence_length=1388,
        matches=1388,
        target_insertions=7661,
        strand=1,
        exons=[
            regions.Exon(start=44506807, stop=44507099),
            regions.Exon(start=44508037, stop=44508309),
            regions.Exon(start=44508727, stop=44508929),
            regions.Exon(start=44509252, stop=44509406),
            regions.Exon(start=44509823, stop=44509927),
            regions.Exon(start=44515497, stop=44515855),
        ]
    ))


@pytest.mark.parametrize('assembly_id,filename,count', [
    ('human', 'data/genome-mapping/results.psl', 275),
])
def test_produces_correct_number_of_results(assembly_id, filename, count):
    assert len(select_hits(assembly_id, filename)) == count


@pytest.mark.parametrize('assembly_id,filename,bad', [
    ('human', 'data/genome-mapping/results.psl', {'URS00000278FE_9606'}),
])
def test_it_does_not_give_data_for_poorly_mapped(assembly_id, filename, bad):
    for result in select_hits(assembly_id, filename):
        assert result.upi not in bad


def test_it_selects_correct_exact_locations():
    val = result_for('human', 'data/genome-mapping/results.psl', 'URS000007D899_9606')
    assert attr.asdict(val) == attr.asdict(gm.Hit(
        assembly_id='human',
        chromosome='21',
        upi='URS000007D899_9606',
        sequence_length=1388,
        matches=1388,
        target_insertions=7661,
        strand=1,
        exons=[
            regions.Exon(start=44506807, stop=44507099),
            regions.Exon(start=44508037, stop=44508309),
            regions.Exon(start=44508727, stop=44508929),
            regions.Exon(start=44509252, stop=44509406),
            regions.Exon(start=44509823, stop=44509927),
            regions.Exon(start=44515497, stop=44515855),
        ]
    ))


def test_selectes_correct_inexact_locations():
    val = results_for('human', 'data/genome-mapping/results.psl', 'URS000093C0D6_9606')
    assert len(val) == 4
    assert attr.asdict(val[0]) == attr.asdict(gm.Hit(
        assembly_id='human',
        chromosome='15',
        upi='URS000093C0D6_9606',
        sequence_length=154,
        matches=153,
        target_insertions=0,
        strand=-1,
        exons=[regions.Exon(start=82478089, stop=82478242)],
    ))
    assert round(val[0].match_fraction, 6) == 0.993506


@pytest.mark.parametrize('assembly_id,filename,upi,region_ids', [
    ('human', 'data/genome-mapping/results.psl', 'URS000007D899_9606', [
        'URS000007D899_9606@21/44506807-44507099,44508037-44508309,44508727-44508929,44509252-44509406,44509823-44509927,44515497-44515855:+'
    ]),
    ('h', 'data/genome-mapping/results.psl', 'URS000093C0D6_9606', [
        'URS000093C0D6_9606@15/82478089-82478242:-',
        'URS000093C0D6_9606@15/84205482-84205576,84205581-84205638:-',
        'URS000093C0D6_9606@15/84467963-84468020,84468025-84468119:+',
        'URS000093C0D6_9606@15/85210081-85210175,85210180-85210237:-',
    ]),
    ('h', 'data/genome-mapping/results.psl', 'URS00002412DD_9606', [
        'URS00002412DD_9606@19/41730635-41730653:+',
        'URS00002412DD_9606@19/42867151-42867169:-',
        'URS00002412DD_9606@19/42944103-42944121:-',
        'URS00002412DD_9606@19/43089440-43089458:-',
        'URS00002412DD_9606@19/43193389-43193407:-',
    ]),
    ('h', 'data/genome-mapping/results.psl', 'URS000073FCEB_9606', [
        'URS000073FCEB_9606@14/20683273-20683309,20683326-20683361:+',
    ]),
    ('h', 'data/genome-mapping/results.psl', 'URS0000584A85_9606', [
        'URS000073FCEB_9606@14/20657464-20657499,20657521-20657557:-',
    ]),
])
def test_can_produce_expected_names(assembly_id, filename, upi, region_ids):
    val = results_for(assembly_id, filename, upi)
    print(val)
    assert [h.name for h in val] == region_ids
