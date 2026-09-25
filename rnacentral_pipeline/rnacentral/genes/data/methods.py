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

import enum

from rnacentral_pipeline.rnacentral.genes.methods.rules import RuleMethod
from rnacentral_pipeline.rnacentral.genes.methods.singletons import SingletonMethod
from rnacentral_pipeline.rnacentral.genes.methods.any_overlap import AnyOverlapMethod


@enum.unique
class Methods(enum.Enum):
    Rules = enum.auto()
    Singleton = enum.auto()
    AnyOverlap = enum.auto()

    @classmethod
    def from_name(cls, raw: str) -> Methods:
        name = raw.lower().replace("_", "").replace("-", "")
        if name == "rules":
            return cls.Rules
        if name == "singleton":
            return cls.Singleton
        if name == "anyoverlap":
            return cls.AnyOverlap
        raise ValueError(f"Unknown method {raw}")

    def handler(self, *args, **kwargs):
        if self is Methods.Rules:
            return RuleMethod(*args, **kwargs)
        if self is Methods.Singleton:
            return SingletonMethod(*args, **kwargs)
        if self is Methods.AnyOverlap:
            return AnyOverlapMethod(*args, **kwargs)
        raise ValueError("Impossible state")
