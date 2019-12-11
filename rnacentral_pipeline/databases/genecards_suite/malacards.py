# -*- coding: utf-8 -*-

"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import typing as ty
import operator as op
import itertools as it

import attr

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.helpers import publications as pub

from .core.data import Context
from .core import parser


def parse(handle, known_handle) -> ty.Iterator[data.Entry]:
    context = Context(
        database='MALACARDS',
        base_url='https://www.malacards.org/card/%s',
        url_data_field='slug',
        gene_field='gene_symbol',
        urs_field='sourceAccession',
        references=[pub.reference(27899610)]
    )
    data = parser.parse(context, handle, known_handle)
    grouped = it.groupby(data, lambda d: d[0].primary_id)
    disease = op.itemgetter('Diseae name')
    for _, items in grouped:
        entries = list(items)
        entry = entries[0][0]
        diseases = sorted({disease(e[1]) for e in entries})
        note_data = entry.note_data
        note_data['diseases'] = diseases
        yield attr.evolve(entry, note_data=note_data)
