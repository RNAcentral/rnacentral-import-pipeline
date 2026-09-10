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


@attr.s(frozen=True)
class Species:
    """
    One piRBase v3.0 species download.

    piRBase names its files by a three letter code and never states a taxid, so
    the mapping lives here. The seventeen already in RNAcentral were taken from
    the taxids the v2 import used; the rest were resolved against NCBI.

    Four needed more than a name lookup. piRBase still calls three of the
    Caenorhabditis by the collection numbers they had before they were
    described, so c26/c31/c32 were resolved through their assemblies
    (JU2190v1, JU2585v1, JU2788v1) to zanzibari, uteleia and sulstoni. psa is
    Plectus sambesii, not the Panagrolaimus the abbreviation suggests.
    """

    code: str = attr.ib(validator=is_a(str))
    name: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int))
    # piRBase is inconsistent about whether the release is in the filename.
    filename: str = attr.ib(validator=is_a(str))
    # Only six species have a curated gold standard set.
    gold: bool = attr.ib(validator=is_a(bool), default=False)


SPECIES = {
    s.code: s
    for s in [
        Species("aca", "Aplysia californica", 6500, "aca.fa.gz"),
        Species("ame", "Ailuropoda melanoleuca", 9646, "ame.v3.0.fa.gz"),
        Species("bgl", "Biomphalaria glabrata", 6526, "bgl.v3.0.fa.gz"),
        Species("bmo", "Bombyx mori", 7091, "bmo.v3.0.fa.gz"),
        Species("bta", "Bos taurus", 9913, "bta.v3.0.fa.gz", True),
        Species("c26", "Caenorhabditis zanzibari", 2306312, "c26.v3.0.fa.gz"),
        Species("c31", "Caenorhabditis uteleia", 2305860, "c31.v3.0.fa.gz"),
        Species("c32", "Caenorhabditis sulstoni", 2305862, "c32.v3.0.fa.gz"),
        Species("cbn", "Caenorhabditis brenneri", 135651, "cbn.v3.0.fa.gz"),
        Species("cbr", "Caenorhabditis briggsae", 6238, "cbr.v3.0.fa.gz"),
        Species("cca", "Caenorhabditis castelli", 1630362, "cca.v3.0.fa.gz"),
        Species("cdo", "Caenorhabditis doughertyi", 1094321, "cdo.v3.0.fa.gz"),
        Species("cel", "Caenorhabditis elegans", 6239, "cel.v3.0.fa.gz"),
        Species("cja", "Callithrix jacchus", 9483, "cja.fa.gz"),
        Species("cma", "Caenorhabditis macrosperma", 1094328, "cma.v3.0.fa.gz"),
        Species("crm", "Caenorhabditis remanei", 31234, "crm.v3.0.fa.gz"),
        Species("cvi", "Caenorhabditis virilis", 1094323, "cvi.v3.0.fa.gz"),
        Species("der", "Drosophila erecta", 7220, "der.fa.gz"),
        Species("dme", "Drosophila melanogaster", 7227, "dme.v3.0.fa.gz", True),
        Species("dpa", "Diploscapter pachys", 2018661, "dpa.v3.0.fa.gz"),
        Species("dre", "Danio rerio", 7955, "dre.fa.gz"),
        Species("dvi", "Drosophila virilis", 7244, "dvi.fa.gz"),
        Species("dya", "Drosophila yakuba", 7245, "dya.fa.gz"),
        Species("eca", "Equus caballus", 9796, "eca.v3.0.fa.gz"),
        Species("gga", "Gallus gallus", 9031, "gga.v3.0.fa.gz"),
        Species("hco", "Haemonchus contortus", 6289, "hco.v3.0.fa.gz"),
        Species("hpo", "Heligmosomoides polygyrus", 6339, "hpo.v3.0.fa.gz"),
        Species("hsa", "Homo sapiens", 9606, "hsa.v3.0.fa.gz", True),
        Species("mfa", "Macaca fascicularis", 9541, "mfa.v3.0.fa.gz", True),
        Species("mml", "Macaca mulatta", 9544, "mml.fa.gz"),
        Species("mmu", "Mus musculus", 10090, "mmu.v3.0.fa.gz", True),
        Species("nbr", "Nippostrongylus brasiliensis", 27835, "nbr.v3.0.fa.gz"),
        Species("nve", "Nematostella vectensis", 45351, "nve.fa.gz"),
        Species("ocu", "Oryctolagus cuniculus", 9986, "ocu.fa.gz"),
        Species("oti", "Oscheius tipulae", 141969, "oti.v3.0.fa.gz"),
        Species("pox", "Poikilolaimus oxycercus", 96659, "pox.v3.0.fa.gz"),
        Species("ppc", "Pristionchus pacificus", 54126, "ppc.v3.0.fa.gz"),
        Species("psa", "Plectus sambesii", 2011161, "psa.v3.0.fa.gz"),
        Species("rno", "Rattus norvegicus", 10116, "rno.fa.gz", True),
        Species("spa", "Scylla paramamosain", 85552, "spa.v3.0.fa.gz"),
        Species("ssc", "Sus scrofa", 9823, "ssc.v3.0.fa.gz"),
        Species("tbe", "Tupaia belangeri", 37347, "tbe.fa.gz"),
        Species("xtr", "Xenopus tropicalis", 8364, "xtr.fa.gz"),
    ]
}


def species(code: str) -> Species:
    if code not in SPECIES:
        raise ValueError(f"Unknown piRBase species: {code}")
    return SPECIES[code]


def urls(base: str) -> ty.Iterable[ty.Tuple[str, str, str]]:
    """
    Every file to fetch, as (code, kind, url). The gold sets are a separate
    download rather than a flag on the full one; they do not share ids with it.
    """

    for s in SPECIES.values():
        yield (s.code, "full", f"{base}/fasta/{s.filename}")
        if s.gold:
            yield (s.code, "gold", f"{base}/fasta/{s.code}.gold.fa.gz")
