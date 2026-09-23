# -*- coding: utf-8 -*-

# pylint: disable=no-member

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

import re
from pathlib import Path

import attr

from rnacentral_pipeline.databases.ensembl.metadata import assemblies as assem

LOAD_SQL = Path("files/import-data/pre-release/000__assemblies.sql")

UCSC = {
    "ucscGenomes": {
        "hg19": {"description": "Feb. 2009 (GRCh37/hg19)"},
        "hg38": {"description": "Dec. 2013 (GRCh38/hg38)"},
    }
}

LINEAGES = {
    9606: "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Chordata; Vertebrata; Mammalia",
    511145: "cellular organisms; Bacteria; Pseudomonadota; Gammaproteobacteria",
    400682: "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Porifera",
    9785: "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Chordata; Vertebrata; Mammalia",
    10228: "cellular organisms; Eukaryota; Opisthokonta; Metazoa; Placozoa",
}


def genome(taxid, accession, name, common_name=None):
    files = {"genes.embl.gz": "e.gz", "genes.gff3.gz": "g.gz"}
    build = {
        "release": "2026_04",
        "paths": {"genebuild": {"files": {"annotations": files}}},
    }
    return {
        "taxid": taxid,
        "common_name": common_name,
        "assemblies": {
            accession: {
                "name": name,
                "level": "chromosome",
                "genebuild_providers": {"ensembl": {"2026_04": build}},
            }
        },
    }


def rows(examples=None, versions=None, **species):
    data = {"species": species}
    return list(
        assem.from_species_json(
            data, frozenset(), examples or {}, UCSC, LINEAGES, versions or {}
        )
    )


def test_builds_a_complete_row_for_a_reference_genome():
    human = genome(9606, "GCA_000001405.29", "GRCh38.p14", common_name="human")
    examples = {"homo_sapiens": {"chromosome": "X", "start": 73819307, "end": 73856333}}
    [row] = rows(examples=examples, Homo_sapiens=human)
    assert attr.asdict(row) == attr.asdict(
        assem.AssemblyInfo(
            assembly_id="GRCh38",
            assembly_full_name="GRCh38.p14",
            gca_accession="GCA_000001405.29",
            assembly_ucsc="hg38",
            common_name="human",
            taxid=9606,
            ensembl_url="homo_sapiens",
            division="EnsemblVertebrates",
            blat_mapping=True,
            example=assem.AssemblyExample(chromosome="X", start=73819307, end=73856333),
        )
    )
    assert row.subdomain == "ensembl.org"


def test_bacteria_get_a_row_now_they_are_imported():
    ecoli = genome(511145, "GCA_000005845.2", "ASM584v2")
    [row] = rows(Escherichia_coli_str_K_12_substr_MG1655=ecoli)
    assert (row.assembly_id, row.division, row.assembly_ucsc) == (
        "ASM584v2",
        "EnsemblBacteria",
        None,
    )


def test_never_writes_the_same_assembly_id_twice():
    """
    Amphimedon and Trichoplax both call their assembly v1.0; assembly_id is the
    table's key, so only the first can have a row.
    """
    got = rows(
        Amphimedon_queenslandica=genome(400682, "GCA_000090795.1", "v1.0"),
        Trichoplax_adhaerens=genome(10228, "GCA_000150275.1", "v1.0"),
    )
    assert [r.ensembl_url for r in got] == ["amphimedon_queenslandica"]


def test_keys_on_the_assembly_id_ensembl_writes_into_its_files():
    """
    species.json names the elephant assembly Loxafr3.0 but the files, and so every
    coordinate the parser writes, say loxAfr3.
    """
    elephant = genome(9785, "GCA_000001905.1", "Loxafr3.0")
    [row] = rows(
        versions={"Loxodonta_africana": "loxAfr3"}, Loxodonta_africana=elephant
    )
    assert (row.assembly_id, row.assembly_full_name) == ("loxAfr3", "Loxafr3.0")


def test_reads_the_genome_version_from_a_gtf_header():
    header = (
        "#!genome-build Loxafr3.0\n#!genome-version loxAfr3\n#!genome-date 2009-07\n"
    )
    assert assem.genome_version(header) == "loxAfr3"
    assert assem.genome_version("1\tensembl\tgene\n") is None


def test_loading_never_blanks_an_existing_ucsc_alias():
    """
    Other databases' hg38/mm10 coordinates are translated through assembly_ucsc,
    so overwriting it with an alias we could not derive would delete them.
    """
    sql = LOAD_SQL.read_text()
    assert re.search(
        r"assembly_ucsc\s*=\s*COALESCE\(\s*EXCLUDED\.assembly_ucsc\s*,\s*ensembl_assembly\.assembly_ucsc\s*\)",
        sql,
        re.IGNORECASE,
    )
