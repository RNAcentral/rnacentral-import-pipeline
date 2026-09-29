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

import re
from pathlib import Path

from rnacentral_pipeline.databases.ensembl import urls

CONFIG = Path(__file__).resolve().parents[3] / "config" / "databases.config"


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


def select(data, legacy=frozenset()):
    return {
        r.species: r
        for r in urls.geneset_urls({"species": data}, base_url="B", legacy=legacy)
    }


def build(release_name, path):
    return release(
        release_name,
        ("genes.embl.gz", f"{path}.embl.gz"),
        ("genes.gff3.gz", f"{path}.gff3.gz"),
    )


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
    assert got.writeable() == (
        "Homo_sapiens",
        "9606",
        "B/p/genes.embl.gz",
        "B/p/genes.gff3.gz",
    )


def test_prefers_ensembl_provider_over_others():
    data = {
        "Sp": species(
            GCA_1=assembly(
                "chromosome",
                {
                    "community": release("2025_01", *STANDARD_FILES),
                    "ensembl": release(
                        "2020_01",
                        ("genes.embl.gz", "e/x.embl.gz"),
                        ("genes.gff3.gz", "e/x.gff3.gz"),
                    ),
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
            GCA_1=assembly(
                "chromosome",
                {
                    "refseq": release(
                        "2026_01",
                        ("genes.embl.gz", "r.embl.gz"),
                        ("genes.gff3.gz", "r.gff3.gz"),
                    )
                },
            ),
            GCA_2=assembly(
                "scaffold",
                {
                    "community": release(
                        "2020_01",
                        ("genes.embl.gz", "c.embl.gz"),
                        ("genes.gff3.gz", "c.gff3.gz"),
                    )
                },
            ),
        )
    }
    got = select(data)["Sp"]
    assert got.embl_url == "B/c.embl.gz"


def test_picks_best_assembly_level_first_within_a_tier():
    data = {
        "Sp": species(
            GCA_1=assembly(
                "scaffold",
                {
                    "ensembl": release(
                        "2026_01",
                        ("genes.embl.gz", "scaf.embl.gz"),
                        ("genes.gff3.gz", "scaf.gff3.gz"),
                    )
                },
            ),
            GCA_2=assembly(
                "chromosome",
                {
                    "ensembl": release(
                        "2020_01",
                        ("genes.embl.gz", "chr.embl.gz"),
                        ("genes.gff3.gz", "chr.gff3.gz"),
                    )
                },
            ),
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
                        **release(
                            "2019_01",
                            ("genes.embl.gz", "old.embl.gz"),
                            ("genes.gff3.gz", "old.gff3.gz"),
                        ),
                        **release(
                            "2026_04",
                            ("genes.embl.gz", "new.embl.gz"),
                            ("genes.gff3.gz", "new.gff3.gz"),
                        ),
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
            GCA_1=assembly(
                "chromosome",
                {"ensembl": release("2026_01", ("genes.embl.gz", "only.embl.gz"))},
            )
        )
    }
    assert "Sp" not in select(data)


def test_reads_the_species_json_ensembl_publishes():
    """
    species.new_ftp_structure.json was retired in August 2026 and now 404s. The
    config value is what the pipeline passes, so it must not drift from the default.
    """
    configured = re.search(r"species_json_url\s*=\s*'([^']+)'", CONFIG.read_text())[1]
    assert configured == urls.DEFAULT_JSON_URL == f"{urls.BASE_URL}/species.json"


def test_keeps_the_assembly_rnacentral_already_uses():
    """
    species.json has no reference flag, and comparing version suffixes across
    accessions picked GRCg6a over GRCg7b, the chicken assembly Ensembl 116 used.
    """
    chicken = {
        "GCA_000002315.5": assembly(
            "chromosome", {"ensembl": build("2022_01", "grcg6a")}
        ),
        "GCA_016699485.1": assembly(
            "chromosome", {"ensembl": build("2022_01", "grcg7b")}
        ),
        "GCA_027557775.1": assembly(
            "chromosome", {"ensembl": build("2023_06", "galgal4")}
        ),
    }
    got = select({"Gallus_gallus": species(**chicken)}, legacy={"GCA_016699485"})
    assert got["Gallus_gallus"].embl_url == "B/grcg7b.embl.gz"


def test_the_assembly_already_used_beats_a_preferred_provider():
    assemblies = {
        "GCA_000000001.1": assembly(
            "chromosome", {"community": build("2018_01", "used")}
        ),
        "GCA_000000002.1": assembly(
            "chromosome", {"ensembl": build("2026_01", "other")}
        ),
    }
    got = select({"Sp": species(**assemblies)}, legacy={"GCA_000000001"})
    assert got["Sp"].embl_url == "B/used.embl.gz"


def test_a_new_species_takes_the_newest_genebuild_not_the_highest_version():
    assemblies = {
        "GCA_000002315.5": assembly("chromosome", {"ensembl": build("2022_01", "old")}),
        "GCA_027557775.1": assembly("chromosome", {"ensembl": build("2023_06", "new")}),
    }
    assert select({"Sp": species(**assemblies)})["Sp"].embl_url == "B/new.embl.gz"


def test_the_used_assembly_list_holds_only_reference_genomes():
    """
    Ensembl 116 also served alternates such as GRCg6a for chicken; listing those
    made them tie with the reference and the first one seen won.
    """
    legacy = urls.legacy_assemblies()
    assert {"GCA_000001405", "GCA_016699485", "GCA_000001735"} <= legacy
    assert "GCA_000002315" not in legacy
    assert all(re.fullmatch(r"GC[AF]_\d{9}", a) for a in legacy)


def test_mrna_rows_carry_the_homology_file_when_there_is_one():
    with_homology = release("2026_04", *STANDARD_FILES)
    with_homology["2026_04"]["paths"]["homologies"] = {
        "files": {"homology_data": {"homology.tsv.gz": "p/homology.tsv.gz"}}
    }
    data = {
        "Homo_sapiens": species(
            **{"GCA_1.1": assembly("chromosome", {"ensembl": with_homology})}
        ),
        "Pan_paniscus": species(
            **{"GCA_2.1": assembly("chromosome", {"ensembl": build("2026_04", "q")})}
        ),
    }
    got = select(data)
    assert got["Homo_sapiens"].writeable(kind="mrna")[-1] == "B/p/homology.tsv.gz"
    assert got["Pan_paniscus"].writeable(kind="mrna")[-1] == ""
    assert len(got["Homo_sapiens"].writeable()) == 4
