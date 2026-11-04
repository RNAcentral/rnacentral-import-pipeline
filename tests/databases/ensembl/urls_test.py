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
        "GCA1v1",
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
    assert all(re.fullmatch(r"GC[AF]_\d{9}(\.\d+)?", a) for a in legacy)
    assert {"GCA_000002765.3", "GCA_016772045.2"} <= legacy


# Genome mapping downloads this file, so it must come from the selected assembly.
def test_carries_the_selected_assemblys_unmasked_genome():
    chosen = assembly("chromosome", {"ensembl": build("2026_04", "new")})
    chosen["assembly"] = {
        "files": {"genome_sequences": {"unmasked.fa.bgz": "g/unmasked.fa.bgz"}}
    }
    other = assembly("scaffold", {"ensembl": build("2026_04", "old")})
    other["assembly"] = {
        "files": {"genome_sequences": {"unmasked.fa.bgz": "o/unmasked.fa.bgz"}}
    }
    got = select({"Danio_rerio": species(**{"GCA_1.1": chosen, "GCA_2.1": other})})
    assert got["Danio_rerio"].genome_url == "B/g/unmasked.fa.bgz"


def test_genome_is_absent_when_species_json_lists_none():
    got = select(
        {
            "Danio_rerio": species(
                **{
                    "GCA_1.1": assembly(
                        "chromosome", {"ensembl": build("2026_04", "x")}
                    )
                }
            )
        }
    )
    assert got["Danio_rerio"].genome_url is None


def test_keeps_the_exact_version_rnacentral_uses():
    assemblies = {
        "GCA_000002765.2": assembly(
            "chromosome", {"community": build("2017_10", "asm276v2")}
        ),
        "GCA_000002765.3": assembly(
            "primary_assembly", {"community": build("2025_08", "with_utrs")}
        ),
    }
    got = select(
        {"Plasmodium_falciparum_3D7": species(**assemblies)},
        legacy={"GCA_000002765", "GCA_000002765.3"},
    )
    assert got["Plasmodium_falciparum_3D7"].embl_url == "B/with_utrs.embl.gz"


def test_one_geneset_per_taxid_when_species_json_lists_a_species_twice():
    """
    species.json lists Aspergillus_nidulans_FGSC_A4 (refseq) and
    Aspergillus_nidulans_fgsc_a4 (community) separately, both taxid 227321.
    """
    data = {
        "Aspergillus_nidulans_FGSC_A4": species(
            GCA_1=assembly("chromosome", {"refseq": build("2015_02", "refseq")})
        ),
        "Aspergillus_nidulans_fgsc_a4": species(
            GCA_2=assembly("chromosome", {"community": build("2015_03", "community")})
        ),
    }
    got = select(data)
    assert list(got) == ["Aspergillus_nidulans_fgsc_a4"]
    assert got["Aspergillus_nidulans_fgsc_a4"].embl_url == "B/community.embl.gz"


def test_one_geneset_per_name_when_species_json_lists_a_strain_under_two_taxids():
    """
    Fusarium_graminearum_PH_1 (taxid 5518) and Fusarium_graminearum_ph_1 (229533)
    are one strain; their names would clash as the assembly's ensembl_url.
    """
    used = {
        "taxid": 5518,
        "assemblies": {
            "GCA_1.1": assembly("chromosome", {"community": build("2025_08", "used")})
        },
    }
    other = {
        "taxid": 229533,
        "assemblies": {
            "GCA_2.3": assembly("chromosome", {"community": build("2013_11", "other")})
        },
    }
    got = select(
        {"Fusarium_graminearum_ph_1": other, "Fusarium_graminearum_PH_1": used},
        legacy={"GCA_1.1"},
    )
    assert list(got) == ["Fusarium_graminearum_PH_1"]


def named(taxid, accession, name):
    entry = assembly("chromosome", {"ensembl": build("2026_04", accession)})
    entry["name"] = name
    return {"taxid": taxid, "assemblies": {accession: entry}}


def ids(data, known=()):
    selected = urls.geneset_urls({"species": data}, base_url="B", known=known)
    return {s.species: s.assembly_id for s in selected}


def test_keeps_the_assembly_id_rnacentral_holds_for_that_gca():
    data = {
        "Plasmodium_falciparum_3D7": named(36329, "GCA_000002765.3", "GCA_000002765")
    }
    got = ids(data, known=[("GCA_000002765.3", "GCA000002765v3")])
    assert got == {"Plasmodium_falciparum_3D7": "GCA000002765v3"}


def test_a_refseq_copy_keeps_the_name_rnacentral_holds_for_its_genbank_twin():
    """
    RefSeq versions count from 1, so GCF_030179895.1 is the copy of GCA_030179895.2.
    """
    data = {
        "Apis_mellifera": named(7460, "GCF_003254395.2", "Amel_HAv3.1"),
        "Drosophila_kikkawai": named(30033, "GCF_030179895.1", "DkikHiC1v2"),
    }
    known = [("GCA_003254395.2", "Amel_HAv3.1"), ("GCA_030179895.2", "DkikHiC1v2")]
    assert ids(data, known) == {
        "Apis_mellifera": "Amel_HAv3.1",
        "Drosophila_kikkawai": "DkikHiC1v2",
    }


def test_a_name_another_genome_holds_gets_a_gca_based_id():
    data = {"Acropora_florida": named(70782, "GCA_014634605.1", "Aflo_1.0")}
    got = ids(data, known=[("GCA_000184785.2", "Aflo_1.0")])
    assert got == {"Acropora_florida": "GCA014634605v1"}


def test_new_genomes_sharing_a_name_get_gca_based_ids():
    data = {
        "Amphimedon_queenslandica": named(400682, "GCA_000090795.1", "v1.0"),
        "Trichoplax_adhaerens": named(10228, "GCA_000150275.1", "v1.0"),
        "Escherichia_coli": named(511145, "GCA_000005845.2", "ASM584v2"),
    }
    assert ids(data) == {
        "Amphimedon_queenslandica": "GCA000090795v1",
        "Trichoplax_adhaerens": "GCA000150275v1",
        "Escherichia_coli": "ASM584v2",
    }


def test_a_new_genome_drops_the_patch_from_its_name():
    data = {"Homo_sapiens": named(9606, "GCA_000001405.29", "GRCh38.p14")}
    assert ids(data) == {"Homo_sapiens": "GRCh38"}


def test_the_assembly_id_list_covers_the_genomes_rnacentral_selects():
    known = dict(urls.known_assemblies())
    assert known["GCA_000001405.29"] == "GRCh38"
    assert known["GCA_000002765.3"] == "GCA000002765v3"
