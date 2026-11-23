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

import collections
import json
import logging
import re
import typing as ty
from pathlib import Path
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

# Assemblies served by Ensembl 116 and Ensembl Genomes 63, the last numbered
# releases. species.json has no reference flag, so a species RNAcentral already
# imports keeps its assembly rather than jumping to a newer or alternate one.
LEGACY_ASSEMBLIES = (
    Path(__file__).resolve().parents[3]
    / "files"
    / "import-data"
    / "ensembl"
    / "legacy-assemblies.txt"
)
# The exact versions RNAcentral selects today, so a species keeps its version too
CURRENT_ASSEMBLIES = LEGACY_ASSEMBLIES.with_name("current-assemblies.txt")
# gca<TAB>assembly_id for every ensembl_assembly row, selected rows first; regenerate
# from production before a run, as Ensembl's own names collide across genomes
ASSEMBLY_IDS = LEGACY_ASSEMBLIES.with_name("assembly-ids.tsv")


@attr.s()
class GenesetUrls:
    species: str = attr.ib(validator=is_a(str))
    taxid: int = attr.ib(validator=is_a(int))
    embl_url: str = attr.ib(validator=is_a(str))
    gff_url: str = attr.ib(validator=is_a(str))
    accession: str = attr.ib(validator=is_a(str))
    homology_url: str = attr.ib(validator=is_a(str), default="")
    genome_url: ty.Optional[str] = attr.ib(default=None)
    assembly_id: ty.Optional[str] = attr.ib(default=None)

    def writeable(self, kind=None) -> ty.Tuple[str, ...]:
        row = (
            self.species,
            str(self.taxid),
            self.embl_url,
            self.gff_url,
            self.assembly_id or "",
        )
        if kind == "mrna":
            return row + (self.homology_url,)
        return row


def _load(location: str):
    """Load the species JSON from an http(s) URL or a local path."""
    if location.startswith("http://") or location.startswith("https://"):
        with urlopen(location) as handle:
            return json.load(handle)
    with open(location, "r") as handle:
        return json.load(handle)


def legacy_assemblies(
    path: Path = LEGACY_ASSEMBLIES, current: Path = CURRENT_ASSEMBLIES
) -> ty.FrozenSet[str]:
    return frozenset(path.read_text().split() + current.read_text().split())


def known_assemblies(path: Path = ASSEMBLY_IDS) -> ty.List[ty.Tuple[str, str]]:
    with path.open() as handle:
        return [tuple(line.rstrip("\n").split("\t")) for line in handle if line.strip()]


def gca_assembly_id(accession: str) -> str:
    # The form RNAcentral already uses for genomes without a usable name
    return accession.replace("_", "").replace(".", "v")


def _assign_assembly_ids(
    data: dict, selected: ty.List[GenesetUrls], known: ty.Iterable[ty.Tuple[str, str]]
) -> None:
    by_gca: ty.Dict[str, str] = {}
    gcas_by_id = collections.defaultdict(set)
    for gca, assembly_id in known:
        by_gca.setdefault(gca, assembly_id)
        gcas_by_id[assembly_id].add(gca)

    number = lambda accession: accession.split("_", 1)[-1].split(".")[0]
    names = {}
    for s in selected:
        name = data["species"][s.species]["assemblies"][s.accession].get("name", "")
        names[s.accession] = re.sub(r"\.p\d+$", "", name)
    shared = collections.Counter(n for a, n in names.items() if a not in by_gca)

    for s in selected:
        name = names[s.accession]
        held = gcas_by_id.get(name, set())
        if s.accession in by_gca:
            s.assembly_id = by_gca[s.accession]
        elif any(number(gca) == number(s.accession) for gca in held):
            # A RefSeq copy of a GenBank assembly RNAcentral holds under this name
            s.assembly_id = name
        elif held or shared[name] > 1 or not name:
            s.assembly_id = gca_assembly_id(s.accession)
        else:
            s.assembly_id = name


def _provider_tier(provider: str) -> int:
    if provider == PREFERRED_PROVIDER:
        return 2
    if provider in NCBI_PROVIDERS:
        return 0
    return 1


def _candidate_key(
    is_current: bool, is_legacy: bool, provider: str, level: str, latest: str
) -> ty.Tuple[bool, bool, int, int, str]:
    # The assembly RNAcentral already uses wins, then provider so a species keeps
    # its Ensembl (or other non-NCBI) build, then assembly level, then newest build.
    return (
        is_current,
        is_legacy,
        _provider_tier(provider),
        LEVEL_RANK.get(level, 0),
        latest,
    )


def _select_geneset(
    species: str, info: dict, base_url: str, legacy: ty.AbstractSet[str]
) -> ty.Optional[ty.Tuple[tuple, GenesetUrls]]:
    """
    Pick the single geneset to import for one species: the most preferred
    (provider, assembly) pair, and that provider's latest release on that
    assembly.
    """

    best = None
    for accession, assembly in info.get("assemblies", {}).items():
        is_legacy = accession.split(".")[0] in legacy
        for provider, releases in assembly.get("genebuild_providers", {}).items():
            key = _candidate_key(
                accession in legacy,
                is_legacy,
                provider,
                assembly.get("level"),
                max(releases),
            )
            if best is None or key > best[0]:
                best = (key, accession, releases)

    if best is None:
        return None

    key, accession, releases = best
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
    homology = (
        release.get("paths", {})
        .get("homologies", {})
        .get("files", {})
        .get("homology_data", {})
        .get("homology.tsv.gz")
    )
    if not embl or not gff:
        LOGGER.warning("Missing EMBL/GFF3 for %s, skipping", species)
        return None

    genome = (
        info["assemblies"][accession]
        .get("assembly", {})
        .get("files", {})
        .get("genome_sequences", {})
        .get("unmasked.fa.bgz")
    )

    return key, GenesetUrls(
        species=species,
        taxid=info["taxid"],
        embl_url=f"{base_url}/{embl}",
        gff_url=f"{base_url}/{gff}",
        accession=accession,
        homology_url=f"{base_url}/{homology}" if homology else "",
        genome_url=f"{base_url}/{genome}" if genome else None,
    )


def _best_by(
    selections: ty.List[ty.Tuple[tuple, GenesetUrls]],
    key: ty.Callable[[GenesetUrls], ty.Hashable],
) -> ty.List[ty.Tuple[tuple, GenesetUrls]]:
    best: ty.Dict[ty.Hashable, ty.Tuple[tuple, GenesetUrls]] = {}
    for selected in selections:
        current = best.get(key(selected[1]))
        if current is None or selected[0] > current[0]:
            best[key(selected[1])] = selected
    return list(best.values())


def geneset_urls(
    data: dict,
    base_url: str = BASE_URL,
    legacy: ty.AbstractSet[str] = frozenset(),
    known: ty.Iterable[ty.Tuple[str, str]] = (),
) -> ty.Iterable[GenesetUrls]:
    selections = []
    for species, info in data["species"].items():
        selected = _select_geneset(species, info, base_url, legacy)
        if selected is not None:
            selections.append(selected)
    # A few species are listed twice, under one taxid or one name in another case.
    selections = _best_by(selections, lambda s: s.taxid)
    chosen = [s for _, s in _best_by(selections, lambda s: s.species.lower())]
    _assign_assembly_ids(data, chosen, known)
    yield from chosen


def urls_for(
    location: str = DEFAULT_JSON_URL, base_url: str = BASE_URL
) -> ty.Iterable[GenesetUrls]:
    """
    Read the unified Ensembl species JSON and yield the geneset EMBL/GFF3 URLs
    to import, one per species.
    """
    yield from geneset_urls(
        _load(location),
        base_url=base_url,
        legacy=legacy_assemblies(),
        known=known_assemblies(),
    )
