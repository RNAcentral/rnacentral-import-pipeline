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

from rnacentral_pipeline.databases import data

from .core.data import Context
from .core import parser


def parse(handle, known_handle) -> ty.Iterator[data.Entry]:
    context = Context(
        database='GENECARDS',
        base_url='https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s',
        url_data_field='GeneCardsSymbol',
        gene_field='GeneCardsSymbol',
        urs_field='URSid',
    )
    yield from parser.parse(context, handle, known_handle)
