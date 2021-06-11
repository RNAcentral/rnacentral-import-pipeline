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

import typing as ty

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.ensembl.vertebrates import urls
from rnacentral_pipeline.databases.ensembl.vertebrates import parser

from rnacentral_pipeline.databases.ensembl.data import FtpInfo


def urls_for(base: str) -> ty.Iterable[FtpInfo]:
    return urls.urls_for(base)


def parse(raw, gff_file, family_file=None) -> ty.Iterable[Entry]:
    return parser.parse(raw, gff_file, family_file=family_file)
