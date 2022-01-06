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

import typing as ty
from pathlib import Path
import operator as op
import itertools as it

import attr

from rnacentral_pipeline.databases.psi_mi import tab
from rnacentral_pipeline.databases.data import Entry, Interaction, InteractionIdentifier

from . import lookup
from . import helpers


def set_interaction_id(interaction: Interaction, index: int) -> Interaction:
    ids = list(interaction.ids)
    value = f"{interaction.urs_taxid}-{index}"
    ids.append(InteractionIdentifier(key="psicquic", value=value, name=""))
    ids = tuple(ids)
    return attr.evolve(interaction, ids=ids)


def parse_interactions(handle: ty.IO) -> ty.Iterable[Interaction]:
    data = tab.parse(handle)
    data = filter(op.methodcaller("involves_rnacentral"), data)
    return data


def parse(path: Path, db_url: str) -> ty.Iterable[Entry]:
    with path.open("r") as raw:
        key = op.attrgetter("urs_taxid")
        interactions = sorted(parse_interactions(raw), key=key)
        mapping = lookup.mapping(db_url, interactions)
        interactions = sorted(interactions, key=key)
        grouped = it.groupby(interactions, key)
        for urs_taxid, current in grouped:
            if urs_taxid not in mapping:
                raise ValueError("Found no sequence info for %s" % urs_taxid)

            info = mapping[urs_taxid]
            current = sorted(current)
            current = [set_interaction_id(inter, i) for i, inter in enumerate(current)]
            entry = helpers.as_entry(urs_taxid, current, info)
            yield entry
