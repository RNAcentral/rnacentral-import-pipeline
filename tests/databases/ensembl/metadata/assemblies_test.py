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

import json
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


def rows(examples=None, known=(), **species):
    data = {"species": species}
    return list(
        assem.from_species_json(
            data, frozenset(), known, examples or {}, UCSC, LINEAGES
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


def test_genomes_sharing_an_assembly_name_each_get_a_row():
    """
    Amphimedon and Trichoplax both call their assembly v1.0; assembly_id is the
    table's key, so each gets an id built from its GCA.
    """
    got = rows(
        Amphimedon_queenslandica=genome(400682, "GCA_000090795.1", "v1.0"),
        Trichoplax_adhaerens=genome(10228, "GCA_000150275.1", "v1.0"),
    )
    assert [r.assembly_id for r in got] == ["GCA000090795v1", "GCA000150275v1"]


def test_names_a_species_by_its_binomial_unless_another_species_shares_it():
    """
    Ensembl called these leishmania_major and escherichia_coli_<strain>; genome
    mapping joins on the name, so two selected genomes must never share one.
    """
    got = rows(
        Leishmania_major_strain_Friedlin=genome(347515, "GCA_000002725.2", "ASM272v2"),
        Escherichia_coli_str_K_12_substr_MG1655=genome(
            511145, "GCA_000005845.2", "ASM584v2"
        ),
        Escherichia_coli_o157_h7=genome(83334, "GCA_004100305.1", "ASM410030v1"),
    )
    assert [r.ensembl_url for r in got] == [
        "leishmania_major",
        "escherichia_coli_str_k_12_substr_mg1655",
        "escherichia_coli_o157_h7",
    ]


def test_plasmodium_falciparum_keeps_blat_mapping_and_its_example():
    """
    Other P. falciparum strains keep 3D7 on its full name, which the BLAT list and
    example locations must use or the existing row loses both.
    """
    examples = json.loads(
        Path("files/import-data/ensembl/example-locations.json").read_text()
    )
    got = rows(
        examples=examples,
        Plasmodium_falciparum_3D7=genome(36329, "GCA_000002765.3", "GCA_000002765"),
        Plasmodium_falciparum_hb3=genome(137071, "GCA_000149665.2", "ASM14966v2"),
    )
    assert (got[0].ensembl_url, got[0].blat_mapping) == (
        "plasmodium_falciparum_3d7",
        True,
    )
    assert got[0].example is not None


def test_anopheles_pest_keeps_blat_mapping_and_its_example():
    """
    species.json also lists a 2022 Anopheles_gambiae build, which takes the plain
    name, so PEST (AgamP4, the assembly RNAcentral uses) keeps its full one.
    """
    examples = json.loads(
        Path("files/import-data/ensembl/example-locations.json").read_text()
    )
    got = rows(
        examples=examples,
        Anopheles_gambiae_str_PEST=genome(180454, "GCA_000005575.1", "AgamP4"),
        Anopheles_gambiae=genome(7165, "GCA_943734735.1", "idAnoGambNW_F1_1"),
    )
    pest = [r for r in got if r.assembly_id == "AgamP4"][0]
    assert (pest.ensembl_url, pest.blat_mapping) == ("anopheles_gambiae_str_pest", True)
    assert pest.example is not None


def test_keeps_the_assembly_id_rnacentral_holds_for_the_genome():
    """
    species.json names the elephant assembly Loxafr3.0; RNAcentral holds it, and
    every coordinate on it, as loxAfr3.
    """
    elephant = genome(9785, "GCA_000001905.1", "Loxafr3.0")
    [row] = rows(known=[("GCA_000001905.1", "loxAfr3")], Loxodonta_africana=elephant)
    assert (row.assembly_id, row.assembly_full_name) == ("loxAfr3", "Loxafr3.0")


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


def test_loading_keeps_an_existing_ensembl_url():
    """
    species.json keys renamed 206 assemblies RNAcentral already holds, e.g.
    leishmania_major became leishmania_major_strain_friedlin.
    """
    sql = LOAD_SQL.read_text()
    assert not re.search(r"ensembl_url\s*=\s*EXCLUDED\.ensembl_url", sql, re.IGNORECASE)


# A missing comma once fused these two into one string, so neither got blat_mapping.
def test_rat_and_yeast_are_blat_genomes():
    assert {"rattus_norvegicus", "saccharomyces_cerevisiae"} <= assem.BLAT_GENOMES
