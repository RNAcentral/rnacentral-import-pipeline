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

import os
import itertools as it

from rnacentral_pipeline.rnacentral.ftp_export import gff3

from tests.helpers import run_with_replacements


def export_text(rna_id, assembly):
    taxid = int(rna_id.split('_')[1])
    path = os.path.join('files', 'ftp-export', 'genome_coordinates', 'query.sql')
    data = run_with_replacements(
        path,
        (':assembly', "'%s'" % assembly),
        (':taxid', "%i" % taxid),
    )
    data = it.ifilter(lambda d: d['rnacentral_id'] == rna_id, data)
    return next(gff3.as_gff3(data))


def test_can_produce_features():
    assert export_text('URS000082BE64_9606', 'hg38') == [
        '3	RNAcentral	transcript	32676136	32676226	.	+	.	ID="LN847466.1:1..91:ncRNA";Name="URS000082BE64";type="snoRNA"',
        '3	RNAcentral	noncoding_exon	32676136	32676226	.	+	.	ID="LN847466.1:1..91:ncRNA_exon1";Name="URS000082BE64";Parent="LN847466.1:1..91:ncRNA";type="snoRNA"',
    ]
