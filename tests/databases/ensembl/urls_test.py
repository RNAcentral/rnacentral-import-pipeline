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

from rnacentral_pipeline.databases.ensembl import urls


def geneset(*files):
    return {"paths": {"genebuild": {"files": {"annotations": dict(files)}}}}


def release(name, *files):
    return {name: {"release": name, **geneset(*files)}}


def assembly(level, providers):
    return {"level": level, "genebuild_providers": providers}


STANDARD_FILES = (
    ("genes.embl.gz", "p/genes.embl.gz"),
    ("genes.gff3.gz", "p/genes.gff3.gz"),
)


def species(**assemblies):
    return {"taxid": 42, "assemblies": assemblies}


def select(data):
    return {r.species: r for r in urls.geneset_urls({"species": data}, base_url="B")}


def test_builds_absolute_urls_from_relative_paths():
    data = {
        "Homo_sapiens": {
            "taxid": 9606,
            "assemblies": {
                "GCA_1.1": assembly(
                    "chromosome",
                    {"ensembl": release("2026_04", *STANDARD_FILES)},
                ),
            },
        }
    }
    got = select(data)["Homo_sapiens"]
    assert got.taxid == 9606
    assert got.embl_url == "B/p/genes.embl.gz"
    assert got.gff_url == "B/p/genes.gff3.gz"
    assert got.writeable() == ("Homo_sapiens", "9606", "B/p/genes.embl.gz", "B/p/genes.gff3.gz")


def test_prefers_ensembl_provider_over_others():
    data = {
        "Sp": species(
            GCA_1=assembly(
                "chromosome",
                {
                    "community": release("2025_01", *STANDARD_FILES),
                    "ensembl": release("2020_01", ("genes.embl.gz", "e/x.embl.gz"), ("genes.gff3.gz", "e/x.gff3.gz")),
                },
            )
        )
    }
    got = select(data)["Sp"]
    assert got.embl_url == "B/e/x.embl.gz"


def test_falls_back_to_ncbi_only_when_nothing_else_and_never_drops():
    data = {
        "Only_refseq": species(
            GCA_1=assembly("scaffold", {"refseq": release("2024_01", *STANDARD_FILES)})
        )
    }
    got = select(data)
    assert "Only_refseq" in got  # no species dropped


def test_prefers_non_ncbi_over_ncbi_even_on_worse_assembly():
    data = {
        "Sp": species(
            GCA_1=assembly("chromosome", {"refseq": release("2026_01", ("genes.embl.gz", "r.embl.gz"), ("genes.gff3.gz", "r.gff3.gz"))}),
            GCA_2=assembly("scaffold", {"community": release("2020_01", ("genes.embl.gz", "c.embl.gz"), ("genes.gff3.gz", "c.gff3.gz"))}),
        )
    }
    got = select(data)["Sp"]
    assert got.embl_url == "B/c.embl.gz"


def test_picks_best_assembly_level_then_version_within_a_tier():
    data = {
        "Sp": species(
            GCA_1=assembly("scaffold", {"ensembl": release("2026_01", ("genes.embl.gz", "scaf.embl.gz"), ("genes.gff3.gz", "scaf.gff3.gz"))}),
            GCA_2=assembly("chromosome", {"ensembl": release("2020_01", ("genes.embl.gz", "chr.embl.gz"), ("genes.gff3.gz", "chr.gff3.gz"))}),
        )
    }
    got = select(data)["Sp"]
    assert got.embl_url == "B/chr.embl.gz"


def test_picks_latest_release():
    data = {
        "Sp": species(
            GCA_1=assembly(
                "chromosome",
                {
                    "ensembl": {
                        **release("2019_01", ("genes.embl.gz", "old.embl.gz"), ("genes.gff3.gz", "old.gff3.gz")),
                        **release("2026_04", ("genes.embl.gz", "new.embl.gz"), ("genes.gff3.gz", "new.gff3.gz")),
                    }
                },
            )
        )
    }
    got = select(data)["Sp"]
    assert got.embl_url == "B/new.embl.gz"


def test_skips_geneset_missing_embl_or_gff():
    data = {
        "Sp": species(
            GCA_1=assembly("chromosome", {"ensembl": release("2026_01", ("genes.embl.gz", "only.embl.gz"))})
        )
    }
    assert "Sp" not in select(data)
