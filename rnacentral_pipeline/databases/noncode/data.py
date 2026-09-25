# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

from __future__ import annotations

import typing as ty

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional


@attr.s(frozen=True)
class Species:
    """
    One NONCODE species download.

    NONCODE names its files by common name or by a three letter code and never
    states a taxid, so the mapping lives here. Every taxid was checked against
    NCBI taxonomy; note gorilla is the species 9593, not the subspecies.
    """

    key: str = attr.ib(validator=is_a(str))
    name: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int))
    # Only human and mouse ship genome coordinates in v6.
    assembly: ty.Optional[str] = attr.ib(validator=optional(is_a(str)), default=None)


SPECIES = {
    s.key: s
    for s in [
        Species("human", "Homo sapiens", 9606, "GRCh38"),
        Species("mouse", "Mus musculus", 10090, "GRCm38"),
        Species("cow", "Bos taurus", 9913),
        Species("rat", "Rattus norvegicus", 10116),
        Species("chicken", "Gallus gallus", 9031),
        Species("celegans", "Caenorhabditis elegans", 6239),
        Species("fruitfly", "Drosophila melanogaster", 7227),
        Species("zebrafish", "Danio rerio", 7955),
        Species("yeast", "Saccharomyces cerevisiae", 4932),
        Species("chimp", "Pan troglodytes", 9598),
        Species("gorilla", "Gorilla gorilla", 9593),
        Species("orangutan", "Pongo abelii", 9601),
        Species("rhesus", "Macaca mulatta", 9544),
        Species("opossum", "Monodelphis domestica", 13616),
        Species("platypus", "Ornithorhynchus anatinus", 9258),
        Species("pig", "Sus scrofa", 9823),
        Species("ATH", "Arabidopsis thaliana", 3702),
        Species("BNA", "Brassica napus", 3708),
        Species("BRA", "Brassica rapa", 3711),
        Species("CQU", "Chenopodium quinoa", 63459),
        Species("CRE", "Chlamydomonas reinhardtii", 3055),
        Species("CSA", "Cucumis sativus", 3659),
        Species("GMA", "Glycine max", 3847),
        Species("GRA", "Gossypium raimondii", 29730),
        Species("MAL", "Malus domestica", 3750),
        Species("MES", "Manihot esculenta", 3983),
        Species("MTR", "Medicago truncatula", 3880),
        Species("MAC", "Musa acuminata", 4641),
        Species("ORU", "Oryza rufipogon", 4529),
        Species("OSA", "Oryza sativa", 4530),
        Species("PPA", "Physcomitrium patens", 3218),
        Species("POP", "Populus trichocarpa", 3694),
        Species("SLY", "Solanum lycopersicum", 4081),
        Species("STU", "Solanum tuberosum", 4113),
        Species("TCA", "Theobroma cacao", 3641),
        Species("TPR", "Trifolium pratense", 57577),
        Species("TAE", "Triticum aestivum", 4565),
        Species("VVI", "Vitis vinifera", 29760),
        Species("ZMA", "Zea mays", 4577),
    ]
}


def species(key: str) -> Species:
    if key not in SPECIES:
        raise ValueError(f"Unknown NONCODE species: {key}")
    return SPECIES[key]
