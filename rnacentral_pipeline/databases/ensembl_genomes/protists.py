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

from rnacentral_pipeline.databases.ensembl_genomes import core
from rnacentral_pipeline.databases.helpers import publications as pubs


def parse(handle):
    context = core.data.Context(
        database='ENSEMBL_PROTISTS',
        references=[pubs.reference('doi:10.1093/nar/gkx1011')],
    )

    return core.parser.parse(context, handle)
