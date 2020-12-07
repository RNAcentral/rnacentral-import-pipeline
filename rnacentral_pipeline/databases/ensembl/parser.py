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

from rnacentral_pipeline.databases.ensembl import fungi
from rnacentral_pipeline.databases.ensembl import metazoa
from rnacentral_pipeline.databases.ensembl import plants
from rnacentral_pipeline.databases.ensembl import protists
from rnacentral_pipeline.databases.ensembl import vertebrates

from rnacentral_pipeline.databases.data import Entry
from rnacentral_pipeline.databases.ensembl.data import Division


def parse(division: Division, *args, **kwargs) -> ty.Iterable[Entry]:
    if division == Division.fungi:
        yield from fungi.parse(*args, **kwargs)
    elif division == Division.metazoa:
        yield from metazoa.parse(*args, **kwargs)
    elif division == Division.plants:
        yield from plants.parse(*args, **kwargs)
    elif division == Division.protists:
        yield from protists.parse(*args, **kwargs)
    elif division == Division.vertebrates:
        yield from vertebrates.parse(*args, **kwargs)
    else:
        raise ValueError(f"Unknown division {division}")
