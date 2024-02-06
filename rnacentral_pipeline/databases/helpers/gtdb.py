# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

import logging
import re
import typing as ty
from functools import lru_cache

from rnacentral_pipeline.databases.helpers import phylogeny as phy

LOGGER = logging.getLogger(__name__)


@lru_cache
def get_inferred_species_taxid(inferred_lineage: str) -> ty.Optional[int]:
    species = re.findall(";s__(.*)", inferred_lineage)[0]
    species = re.sub(r"_\w+", "", species)
    if len(species) == 0:
        ## There is no species in the inferred phylogeny, so try next layer up
        return None
    taxid = None

    try:
        taxid = phy.taxid(species)
    except:
        LOGGER.warning(f"{species} doesn't have a taxid, trying something else")
        species_name = species.split()[0]
        env_sample_name = f"uncultured {species_name} sp."
        try:
            taxid = phy.taxid(env_sample_name)
        except:
            LOGGER.warning(
                f"Uncultured species {env_sample_name} doesn't have a taxid, trying something else"
            )

    return taxid


@lru_cache
def get_inferred_genus_taxid(inferred_lineage: str) -> ty.Optional[int]:
    genus = re.findall(";g__(.*);", inferred_lineage)[0]
    if len(genus) == 0:
        ## There is no genus in the inferred phylogeny, so try next layer up
        return None

    genus = genus.split("_")[0]
    try:
        return phy.taxid(f"uncultured {genus} sp.")
    except:
        LOGGER.warning(f"{genus} doesn't have a taxid, trying something else")
    return None


@lru_cache
def get_inferred_family_taxid(inferred_lineage: str) -> ty.Optional[int]:
    family = re.findall(";f__(.*);g", inferred_lineage)[0]
    if len(family) == 0:
        ## There is no family in the inferred phylogeny, so try next layer up
        return None

    try:
        return phy.taxid(f"{family} bacterium")
    except:
        LOGGER.warning(f"{family} doesn't have a taxid, trying something else")
    return None


@lru_cache
def get_inferred_order_taxid(inferred_lineage: str) -> ty.Optional[int]:
    order = re.findall(";o__(.*);f", inferred_lineage)[0]
    if len(order) == 0:
        ## There is no order in the inferred phylogeny, so try next layer up
        return None
    try:
        return phy.taxid(f"{order}")
    except:
        LOGGER.warning(f"{order} doesn't have a taxid, trying something else")
    return None


@lru_cache
def phylogeny_to_taxid(lineage: str) -> int:
    getters = [
        get_inferred_species_taxid,
        get_inferred_genus_taxid,
        get_inferred_family_taxid,
        get_inferred_order_taxid,
    ]
    for getter in getters:
        if value := getter(lineage):
            return value
    raise ValueError(f"Could not get taxid for `${lineage}`")
