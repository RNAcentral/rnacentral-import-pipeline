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

import re

import attr

from rnacentral_pipeline.databases.ensembl.genomes import parser
from rnacentral_pipeline.databases.ensembl.genomes.data import Context
from rnacentral_pipeline.databases.helpers import publications as pubs


def correct_rna_type(entry):
    if re.match(r'^LMJF_\d+_snoRNA_?\d+$', entry.gene):
        return attr.evolve(entry, rna_type='snoRNA')

    if re.match(r'^LMJF_\d+_snRNA_\d+$', entry.gene):
        return attr.evolve(entry, rna_type='snRNA')

    if entry.rna_type == 'snRNA' and 'snoRNA' in entry.description:
        return attr.evolve(entry, rna_type='snoRNA')

    return entry


def parse(handle, gff_file):
    context = Context.build(
        'ENSEMBL_PROTISTS',
        [pubs.reference('doi:10.1093/nar/gkx1011')],
        gff_file,
    )

    return map(correct_rna_type, parser.parse(context, handle))
