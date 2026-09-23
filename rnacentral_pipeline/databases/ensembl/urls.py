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

import json
import logging
import typing as ty
from urllib.request import urlopen

import attr
from attr.validators import instance_of as is_a

LOGGER = logging.getLogger(__name__)

# Base of the unified Ensembl organisms FTP. Every path in the species JSON is
# relative to this.
BASE_URL = "https://ftp.ebi.ac.uk/pub/ensemblorganisms"
DEFAULT_JSON_URL = f"{BASE_URL}/species.json"

# Import the primary annotation for every species (no species is dropped). When
# a species has more than one gene build we prefer Ensembl's own, then any other
# non-NCBI producer, and fall back to the NCBI-derived builds (which RNAcentral
# also imports via NCBI/RefSeq) only when nothing else is available. See the
# update/Ensembl-import decisions.
NCBI_PROVIDERS = frozenset({"refseq", "genbank"})
PREFERRED_PROVIDER = "ensembl"

# Rank assembly quality to pick the single best assembly for a species. Higher
# is better; unknown levels sort last.
LEVEL_RANK = {
    "chromosome_group": 6,
    "chromosome": 5,
    "primary_assembly": 4,
    "supercontig": 3,
    "scaffold": 2,
    "contig": 1,
}


@attr.s()
class GenesetUrls:
    species: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int))
    embl_url: str = attr.ib(validator=is_a(str))
    gff_url: str = attr.ib(validator=is_a(str))

    def writeable(self, kind=None) -> ty.Tuple[str, str, str, str]:
        return (self.species, str(self.taxid), self.embl_url, self.gff_url)


def _load(location: str):
    """Load the species JSON from an http(s) URL or a local path."""
    if location.startswith("http://") or location.startswith("https://"):
        with urlopen(location) as handle:
            return json.load(handle)
    with open(location, "r") as handle:
        return json.load(handle)


def _assembly_version(accession: str) -> int:
    try:
        return int(accession.rsplit(".", 1)[1])
    except (IndexError, ValueError):
        return 0


def _provider_tier(provider: str) -> int:
    if provider == PREFERRED_PROVIDER:
        return 2
    if provider in NCBI_PROVIDERS:
        return 0
    return 1


def _candidate_key(level: str, version: int, provider: str) -> ty.Tuple[int, int, int]:
    # Provider preference wins first so a species keeps its Ensembl (or other
    # non-NCBI) build; then the best assembly (level, then version).
    return (_provider_tier(provider), LEVEL_RANK.get(level, 0), version)


def _select_geneset(
    species: str, info: dict, base_url: str
) -> ty.Optional[GenesetUrls]:
    """
    Pick the single geneset to import for one species: the most preferred
    (provider, assembly) pair, and that provider's latest release on that
    assembly.
    """

    best = None
    for accession, assembly in info.get("assemblies", {}).items():
        version = _assembly_version(accession)
        for provider, releases in assembly.get("genebuild_providers", {}).items():
            key = _candidate_key(assembly.get("level"), version, provider)
            if best is None or key > best[0]:
                best = (key, releases)

    if best is None:
        return None

    releases = best[1]
    # Release keys are zero-padded YYYY_MM, so a string max is the latest.
    _, release = max(releases.items(), key=lambda kv: kv[0])

    annotations = (
        release.get("paths", {})
        .get("genebuild", {})
        .get("files", {})
        .get("annotations", {})
    )
    embl = annotations.get("genes.embl.gz")
    gff = annotations.get("genes.gff3.gz")
    if not embl or not gff:
        LOGGER.warning("Missing EMBL/GFF3 for %s, skipping", species)
        return None

    return GenesetUrls(
        species=species,
        taxid=info["taxid"],
        embl_url=f"{base_url}/{embl}",
        gff_url=f"{base_url}/{gff}",
    )


def geneset_urls(data: dict, base_url: str = BASE_URL) -> ty.Iterable[GenesetUrls]:
    for species, info in data["species"].items():
        selected = _select_geneset(species, info, base_url)
        if selected is not None:
            yield selected


def urls_for(
    location: str = DEFAULT_JSON_URL, base_url: str = BASE_URL
) -> ty.Iterable[GenesetUrls]:
    """
    Read the unified Ensembl species JSON and yield the geneset EMBL/GFF3 URLs
    to import, one per species.
    """
    yield from geneset_urls(_load(location), base_url=base_url)
