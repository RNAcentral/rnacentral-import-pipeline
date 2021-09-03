# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.ensembl import fungi
from rnacentral_pipeline.databases.ensembl import metazoa
from rnacentral_pipeline.databases.ensembl import plants
from rnacentral_pipeline.databases.ensembl import protists
from rnacentral_pipeline.databases.ensembl import vertebrates

from rnacentral_pipeline.databases.ensembl.data import Division, Pseudogene


def parse(division: Division, handle: ty.IO) -> ty.Iterable[Pseudogene]:
    if division == Division.fungi:
        yield from fungi.pseudogenes(handle)
    elif division == Division.metazoa:
        yield from metazoa.pseudogenes(handle)
    elif division == Division.plants:
        yield from plants.pseudogenes(handle)
    elif division == Division.protists:
        yield from protists.pseudogenes(handle)
    elif division == Division.vertebrates:
        yield from vertebrates.pseudogenes(handle)
    else:
        raise ValueError(f"Unknown division {division}")
