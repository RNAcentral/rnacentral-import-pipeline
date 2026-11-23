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

import collections
import json
import logging
import re
import typing as ty
from pathlib import Path

import attr
from attr.validators import instance_of as is_a
from attr.validators import optional

from rnacentral_pipeline import schemas
from rnacentral_pipeline.databases.ensembl import urls
from rnacentral_pipeline.parquet_writers import row_writer

BLAT_GENOMES = {
    "anopheles_gambiae_str_pest",
    "arabidopsis_thaliana",
    "bombyx_mori",
    "caenorhabditis_elegans",
    "dictyostelium_discoideum",
    "drosophila_melanogaster",
    "homo_sapiens",
    "mus_musculus",
    "plasmodium_falciparum_3d7",
    "rattus_norvegicus",
    "saccharomyces_cerevisiae",
    "schizosaccharomyces_pombe",
}

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
    blat_mapping = attr.ib(validator=is_a(bool), converter=bool)
    example = attr.ib(validator=optional(is_a(AssemblyExample)))

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
            chromosome,
            start,
            end,
            int(self.blat_mapping),
        ]


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


def ensembl_urls(species: ty.Iterable[str]) -> ty.Dict[str, str]:
    # Like Ensembl's old names: the binomial, unless strains of it need the full name.
    names = [s.lower() for s in species]
    binomial = {n: "_".join(n.split("_")[:2]) for n in names}
    counts = collections.Counter(binomial.values())
    return {n: binomial[n] if counts[binomial[n]] == 1 else n for n in names}


def from_species_json(
    data: dict,
    legacy: ty.AbstractSet[str],
    known: ty.Iterable[ty.Tuple[str, str]],
    examples: dict,
    ucsc: dict,
) -> ty.Iterable[AssemblyInfo]:
    """
    One row per geneset the import selects, keyed by the assembly id the selection
    gives it, which the parser also puts on every coordinate.
    """
    aliases = ucsc_aliases(ucsc)
    selections = list(urls.geneset_urls(data, legacy=legacy, known=known))
    names = ensembl_urls(s.species for s in selections)
    seen = set()
    for selected in selections:
        info = data["species"][selected.species]
        name = info["assemblies"][selected.accession]["name"]
        assembly_id = selected.assembly_id
        if assembly_id in seen:
            LOGGER.warning(
                "Assembly %s already used, skipping %s", assembly_id, selected.species
            )
            continue
        seen.add(assembly_id)

        url = names[selected.species.lower()]
        example = examples.get(url)
        yield AssemblyInfo(
            assembly_id=assembly_id,
            assembly_full_name=name,
            gca_accession=selected.accession,
            assembly_ucsc=aliases.get(assembly_id),
            common_name=info.get("common_name"),
            taxid=reconcile_taxids(info["taxid"]),
            ensembl_url=url,
            blat_mapping=url in BLAT_GENOMES,
            example=AssemblyExample(
                chromosome=str(example["chromosome"]),
                start=example["start"],
                end=example["end"],
            )
            if example
            else None,
        )


def write(location, example_file, ucsc_file, output):
    """
    Write the ensembl_assembly rows for the genesets selected from species.json.
    """

    data = urls._load(location)
    rows = from_species_json(
        data,
        urls.legacy_assemblies(),
        urls.known_assemblies(),
        json.load(example_file),
        json.load(ucsc_file),
    )
    with row_writer(Path(output), schemas.ASSEMBLIES) as writer:
        writer.writerows(row.writeable() for row in rows)
