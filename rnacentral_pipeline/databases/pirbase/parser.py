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

import ijson
from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.generic import v1


def parse(handle) -> ty.Iterable[data.Entry]:
    for entry in ijson.items(handle, "data"):
        yield v1.as_ncrna(entry)
