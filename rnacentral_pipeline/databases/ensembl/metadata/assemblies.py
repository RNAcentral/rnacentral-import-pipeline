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

import json
import logging
import re
import typing as ty
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import attr
import psycopg2 as pg
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline import schemas
from rnacentral_pipeline.databases.ensembl import urls
from rnacentral_pipeline.parquet_writers import row_writer

BLAT_GENOMES = {
    "anopheles_gambiae",
    "arabidopsis_thaliana",
    "bombyx_mori",
    "caenorhabditis_elegans",
    "dictyostelium_discoideum",
    "drosophila_melanogaster",
    "homo_sapiens",
    "mus_musculus",
    "plasmodium_falciparum",
    "rattus_norvegicus" "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe",
}

# First match wins, so Vertebrata must come before Metazoa. Anything else is a
# protist, as in the FTP export's kingdom split.
DIVISIONS = [
    ("Vertebrata", "EnsemblVertebrates"),
    ("Viridiplantae", "EnsemblPlants"),
    ("Fungi", "EnsemblFungi"),
    ("Metazoa", "EnsemblMetazoa"),
    ("Bacteria", "EnsemblBacteria"),
    ("Archaea", "EnsemblBacteria"),
]

LOGGER = logging.getLogger(__name__)


def reconcile_taxids(taxid):
    """
    Sometimes Ensembl taxid and ENA/Expert Database taxid do not match,
    so to reconcile the differences, taxids in the ensembl_assembly table
    are overriden to match other data.
    """

    taxid = str(taxid)
    if taxid == "284812":  # Ensembl assembly for Schizosaccharomyces pombe
        return 4896  # Pombase and ENA xrefs for Schizosaccharomyces pombe
    return int(taxid)


class InvalidDomain(Exception):
    """
    Raised when we cannot compute a URL for a domain.
    """

    pass


@attr.s()
class AssemblyExample(object):
    chromosome = attr.ib(validator=is_a(str))
    start = attr.ib(validator=is_a(int))
    end = attr.ib(validator=is_a(int))


@attr.s()
class AssemblyInfo(object):
    assembly_id = attr.ib(validator=is_a(str))
    assembly_full_name = attr.ib(validator=is_a(str))
    gca_accession = attr.ib(validator=optional(is_a(str)))
    assembly_ucsc = attr.ib(validator=optional(is_a(str)))
    common_name = attr.ib(validator=optional(is_a(str)))
    taxid = attr.ib(validator=is_a(int))
    ensembl_url = attr.ib(validator=is_a(str))
    division = attr.ib(validator=is_a(str))
    blat_mapping = attr.ib(validator=is_a(bool), converter=bool)
    example = attr.ib(validator=optional(is_a(AssemblyExample)))

    @property
    def subdomain(self):
        """Given E! division, returns E!/E! Genomes url."""

        if self.division == "Ensembl":
            return "ensembl.org"
        if self.division == "EnsemblPlants":
            return "plants.ensembl.org"
        if self.division == "EnsemblMetazoa":
            return "metazoa.ensembl.org"
        if self.division == "EnsemblBacteria":
            return "bacteria.ensembl.org"
        if self.division == "EnsemblFungi":
            return "fungi.ensembl.org"
        if self.division == "EnsemblProtists":
            return "protists.ensembl.org"
        if self.division == "EnsemblVertebrates":
            return "ensembl.org"
        raise InvalidDomain(self.division)

    def writeable(self):
        chromosome = None
        start = None
        end = None
        if self.example:
            chromosome = self.example.chromosome
            start = self.example.start
            end = self.example.end

        return [
            self.assembly_id,
            self.assembly_full_name,
            self.gca_accession,
            self.assembly_ucsc,
            self.common_name,
            self.taxid,
            self.ensembl_url,
            self.division,
            self.subdomain,
            chromosome,
            start,
            end,
            int(self.blat_mapping),
        ]


def division(lineage: str) -> str:
    for marker, name in DIVISIONS:
        if marker in lineage:
            return name
    return "EnsemblProtists"


def ucsc_aliases(ucsc: dict) -> ty.Dict[str, str]:
    """
    Map assembly name to UCSC database, read from descriptions like
    "Dec. 2013 (GRCh38/hg38)". UCSC's GCA accessions are the unpatched
    versions and are shared by GRCh37 and GRCh38, so they cannot be used.
    """
    aliases = {}
    for db, info in ucsc["ucscGenomes"].items():
        found = re.search(r"\(([^()/]+)/" + re.escape(db) + r"\)", info["description"])
        if found:
            aliases.setdefault(found[1], db)
    return aliases


def genome_version(header: str) -> ty.Optional[str]:
    for line in header.splitlines():
        if line.startswith("#!genome-version"):
            return line.split(" ", 1)[1].strip()
    return None


def fetch_genome_version(selected: urls.GenesetUrls) -> ty.Optional[str]:
    """
    The GTF carries the same header as the GFF3 but without the sequence-region
    lines, so the first couple of KB always hold it.
    """
    url = selected.gff_url.replace("genes.gff3.gz", "genes.gtf.gz")
    request = urllib.request.Request(url, headers={"Range": "bytes=0-2047"})
    for _ in range(3):
        try:
            raw = urllib.request.urlopen(request, timeout=60).read()
            text = zlib.decompressobj(16 + zlib.MAX_WBITS).decompress(raw)
            return genome_version(text.decode(errors="replace"))
        except Exception as err:
            LOGGER.warning("Could not read %s: %s", url, err)
    return None


def from_species_json(
    data: dict,
    legacy: ty.AbstractSet[str],
    examples: dict,
    ucsc: dict,
    lineages: ty.Dict[int, str],
    versions: ty.Dict[str, str],
) -> ty.Iterable[AssemblyInfo]:
    """
    One row per geneset the import selects, keyed by the assembly id Ensembl
    writes into the files, which is what the parser puts on every coordinate.
    species.json's name without its patch suffix agrees ~96% of the time, so
    it is only the fallback (elephant: Loxafr3.0 there, loxAfr3 in the files).
    """
    aliases = ucsc_aliases(ucsc)
    seen = set()
    for selected in urls.geneset_urls(data, legacy=legacy):
        info = data["species"][selected.species]
        name = info["assemblies"][selected.accession]["name"]
        assembly_id = versions.get(selected.species)
        if not assembly_id:
            assembly_id = re.sub(r"\.p\d+$", "", name)
            LOGGER.warning(
                "No genome version for %s, using %s", selected.species, assembly_id
            )
        if assembly_id in seen:
            LOGGER.warning(
                "Assembly %s already used, skipping %s", assembly_id, selected.species
            )
            continue
        seen.add(assembly_id)

        url = selected.species.lower()
        example = examples.get(url)
        yield AssemblyInfo(
            assembly_id=assembly_id,
            assembly_full_name=name,
            gca_accession=selected.accession,
            assembly_ucsc=aliases.get(assembly_id),
            common_name=info.get("common_name"),
            taxid=reconcile_taxids(info["taxid"]),
            ensembl_url=url,
            division=division(lineages.get(info["taxid"], "")),
            blat_mapping=url in BLAT_GENOMES,
            example=AssemblyExample(
                chromosome=str(example["chromosome"]),
                start=example["start"],
                end=example["end"],
            )
            if example
            else None,
        )


def load_lineages(db_url: str, taxids: ty.List[int]) -> ty.Dict[int, str]:
    with pg.connect(db_url) as conn, conn.cursor() as cur:
        cur.execute(
            "select id, lineage from rnc_taxonomy where id = any(%s)", (taxids,)
        )
        return dict(cur.fetchall())


def write(location, example_file, ucsc_file, output, db_url):
    """
    Write the ensembl_assembly rows for the genesets selected from species.json.
    """

    data = urls._load(location)
    legacy = urls.legacy_assemblies()
    selected = list(urls.geneset_urls(data, legacy=legacy))
    with ThreadPoolExecutor(10) as pool:
        found = pool.map(fetch_genome_version, selected)
        versions = {s.species: v for s, v in zip(selected, found) if v}

    taxids = [info["taxid"] for info in data["species"].values()]
    rows = from_species_json(
        data,
        legacy,
        json.load(example_file),
        json.load(ucsc_file),
        load_lineages(db_url, taxids),
        versions,
    )
    with row_writer(Path(output), schemas.ASSEMBLIES) as writer:
        writer.writerows(row.writeable() for row in rows)
