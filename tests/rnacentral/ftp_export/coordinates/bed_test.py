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

from rnacentral_pipeline.rnacentral.ftp_export.coordinates import bed

from .helpers import fetch_coord


def fetch_data(rna_id, assembly):
    data = [bed.BedEntry.from_region(c) for c in fetch_coord(rna_id, assembly)]
    assert len(data) == 1
    return data[0]


def test_can_build_bed_from_region():
    data = fetch_data('URS000082BE64_9606', "GRCh38")
    assert attr.asdict(data) == attr.asdict(bed.BedEntry(
        chromosome='3',
        rna_id='URS000082BE64_9606',
        blocks=[bed.BedBlock(start=32676136, stop=32676226)],
        strand=1,
        rna_type='snoRNA',
        databases='snOPY',
    ))


def test_can_build_entry_with_several_databases():
    data = fetch_data('URS0000368518_9606', "GRCh38")
    assert attr.asdict(data) == attr.asdict(bed.BedEntry(
        chromosome='2',
        rna_id='URS0000368518_9606',
        blocks=[bed.BedBlock(start=37208875, stop=37212677)],
        strand=1,
        rna_type='lncRNA',
        databases='Ensembl,GENCODE',
    ))
