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
    species = re.findall(";s__(.*)", inferred_lineage)
    if not species:
        return None
    species = re.sub(r"_\w+", "", species[0])
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
    genus = re.findall(";g__(.*);", inferred_lineage)
    if not genus or not genus[0]:
        ## There is no genus in the inferred phylogeny, so try next layer up
        return None

    genus = genus[0].split("_")[0]
    try:
        return phy.taxid(f"uncultured {genus} sp.")
    except:
        LOGGER.warning(f"{genus} doesn't have a taxid, trying something else")
    return None


@lru_cache
def get_inferred_family_taxid(inferred_lineage: str) -> ty.Optional[int]:
    family = re.findall(r";f__(.*)[;g|$]", inferred_lineage)
    if not family or not family[0]:
        ## There is no family in the inferred phylogeny, so try next layer up
        return None

    try:
        return phy.taxid(f"{family} bacterium")
    except:
        LOGGER.warning(f"{family} doesn't have a taxid, trying something else")
    return None


@lru_cache
def get_inferred_order_taxid(inferred_lineage: str) -> ty.Optional[int]:
    """
    Get the taxid for the order entry, if any in the given lineage. If no order
    entry can be found then None is returned.

    >>> get_inferred_order_taxid("d__Bacteria,p__UBA6262,c__UBA6262")
    None
    >>> get_inferred_order_taxid("Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae")
    178369
    >>> get_inferred_order_taxid("Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida")
    178369
    """
    order = re.findall("[;,]o__(.*)[;,]f", inferred_lineage)
    if not order or not order[0]:
        ## There is no order in the inferred phylogeny, so try next layer up
        return None
    try:
        return phy.taxid(f"{order}")
    except:
        LOGGER.warning(f"{order} doesn't have a taxid, trying something else")
    return None


@lru_cache
def get_by_part(
    lineage: str, bacteria_fallback: str | int = "bacterium"
) -> ty.Optional[int]:
    """
    Tries to fetch the NCBI taxid for a lineage from GTDB. The lineage can be
    ';' or ',' separated and handles the 's__' style level indicators. This will try
    one level at a time, starting with the most specific. In the case this has
    fallen back to 'Bacteria', this will use fetch the value value for the bacteria_fallback,
    name, which defaults to 'bacterium', eg:

    https://www.ncbi.nlm.nih.gov/datasets/taxonomy/1869227/.

    If no taxonomic id can be found, None is returned.

    >>> get_by_part("Chromatophores,d__Eukaryota,p__Cercozoa,c__Imbricatea,o__Euglyphida,f__Paulinellidae")
    39715
    >>> get_by_part("Chromatophores; d__Eukaryota; p__Cercozoa; c__Imbricatea; o__Euglyphida; f__Paulinellidae")
    39715
    >>> get_by_part("d__Bacteria,p__UBA6262,c__UBA6262,o__WVXT01,f__WVXT01,WVXT01__sp009619095")
    1869227
    >>> get_by_part("d__Bacteria,p__UBA6262,c__UBA6262,o__WVXT01,f__WVXT01,WVXT01__sp009619095", bacteria_fallback=2)
    2
    >>> get_by_part("d__Bacteria,p__UBA6262,c__UBA6262,o__WVXT01,f__WVXT01,WVXT01__sp009619095", bacteria_fallback="bacteria")
    2
    >>> get_by_part("p__UBA6262,c__UBA6262,o__WVXT01,f__WVXT01,WVXT01__sp009619095")
    None
    """
    parts = re.split(r"[;,]\s*", lineage)
    for part in reversed(parts):
        name = part
        if match := re.match(r"^(\w)__(.+)", part):
            name = match.group(2)
        if name.lower() == "bacteria" and bacteria_fallback:
            if isinstance(bacteria_fallback, int):
                return bacteria_fallback
            return phy.taxid(bacteria_fallback)
        try:
            return phy.taxid(name.lower())
        except:
            LOGGER.info("Failed to get taxid for %s (%s)", name, part)
    return None


@lru_cache
def phylogeny_to_taxid(lineage: str) -> int:
    """
    Tries to fetch the phylogeny using a variety of strategies. This will try
    the get_by_part function, and then all other strategies in this module.
    """
    getters = [
        get_by_part,
        get_inferred_species_taxid,
        get_inferred_genus_taxid,
        get_inferred_family_taxid,
        get_inferred_order_taxid,
    ]
    for getter in getters:
        if value := getter(lineage):
            return value
    raise ValueError(f"Could not get taxid for `{lineage}`")
