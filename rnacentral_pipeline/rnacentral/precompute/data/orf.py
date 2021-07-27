from __future__ import annotations

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

import attr
from attr.validators import instance_of as is_a


@attr.s()
class OrfInfo:
    sources: ty.List[str] = attr.ib(validator=is_a(list))

    @classmethod
    def build(cls, raw) -> OrfInfo:
        return cls(**raw)

    def all_sources(self) -> ty.List[str]:
        sources = set(self.sources)
        if 'cpat' in sources:
            sources.remove('cpat')
            sources.add('CPAT')
        return sorted(sources)
