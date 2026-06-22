# -*- coding: utf-8 -*-

"""
Copyright [2009-2024] EMBL-European Bioinformatics Institute
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

from io import StringIO

import pytest

from rnacentral_pipeline.databases.ncbi import taxonomy
from rnacentral_pipeline.databases.helpers import phylogeny as phy

# --- nodes.dmp fixtures ---
# Real format: taxid \t|\t parent_taxid \t|\t rank \t|\t ...more fields... \t|\n
NODES_DMP = (
    "9606\t|\t9605\t|\tspecies\t|\n"
    "10090\t|\t10088\t|\tspecies\t|\n"
    "9443\t|\t314146\t|\torder\t|\n"
)

# --- UniProt STATS file fixtures ---
# 15 preamble lines, then tab-separated data with taxid in column 1 (0-indexed)
STATS_HEADER = "".join(f"preamble line {i}\n" for i in range(15))
STATS_DATA = (
    "UP000005640\t9606\tHomo sapiens\t20600\n"
    "UP000000589\t10090\tMus musculus\t17100\n"
)

# --- fullnamelineage.dmp fixture ---
# Format: taxid \t|\t name \t|\t lineage \t|\n
LINEAGE_DMP = (
    "9606\t|\tHomo sapiens\t|\tEukaryota; Metazoa; Mammalia; Primates; \t|\n"
    "9443\t|\tPrimates\t|\tEukaryota; Metazoa; Mammalia; \t|\n"
)

# --- names.dmp fixture ---
# Format: taxid \t|\t name \t|\t unique_name \t|\t name_class \t|\n
NAMES_DMP = (
    "9606\t|\tHomo sapiens\t|\t\t|\tscientific name\t|\n"
    "9606\t|\thuman\t|\t\t|\tgenbank common name\t|\n"
    "9443\t|\tPrimates\t|\t\t|\tscientific name\t|\n"
)

# --- merged.dmp fixture (empty — no merged taxids) ---
MERGED_DMP = ""

# --- delnodes.dmp fixture ---
# Real format: deleted_taxid \t|\n (just a bare id per line). 2 is deleted; 9606
# is included to check that a still-live node is never emitted as deleted.
DELNODES_DMP = "2\t|\n" "9606\t|\n"

# --- delnodes.dmp fixture (empty — no deleted taxids) ---
DELNODES_DMP_EMPTY = ""


def test_parse_nodes():
    result = taxonomy.parse_nodes(StringIO(NODES_DMP))
    assert result == {"9606": "species", "10090": "species", "9443": "order"}


def test_parse_delnodes():
    result = taxonomy.parse_delnodes(StringIO(DELNODES_DMP))
    assert result == {"2", "9606"}


def test_parse_delnodes_empty():
    result = taxonomy.parse_delnodes(StringIO(DELNODES_DMP_EMPTY))
    assert result == set()


def test_parse_ref_proteomes():
    result = taxonomy.parse_ref_proteomes(StringIO(STATS_HEADER + STATS_DATA))
    assert result == {9606, 10090}


def test_parse_ref_proteomes_empty():
    result = taxonomy.parse_ref_proteomes(StringIO(STATS_HEADER))
    assert result == set()


def test_taxonomy_entry_build_with_rank_and_proteome():
    entry = taxonomy.TaxonomyEntry.build(
        entry=["9606", "Homo sapiens", "Eukaryota; Metazoa; "],
        names=[("9606", "Homo sapiens", "", "scientific name")],
        rank="species",
        reference_proteome=True,
    )
    assert entry.rank == "species"
    assert entry.reference_proteome is True
    assert entry.tax_id == 9606


def test_taxonomy_entry_build_defaults():
    entry = taxonomy.TaxonomyEntry.build(
        entry=["9606", "Homo sapiens", "Eukaryota; Metazoa; "],
        names=[("9606", "Homo sapiens", "", "scientific name")],
    )
    assert entry.rank == ""
    assert entry.reference_proteome is False


def test_taxonomy_entry_writeable_includes_rank_and_proteome():
    entry = taxonomy.TaxonomyEntry.build(
        entry=["9606", "Homo sapiens", "Eukaryota; Metazoa; "],
        names=[("9606", "Homo sapiens", "", "scientific name")],
        rank="species",
        reference_proteome=True,
    )
    rows = list(entry.writeable())
    assert len(rows) == 1
    row = rows[0]
    assert row[5] == "species"
    assert row[6] is True
    assert row[7] is False


def test_taxonomy_entry_build_deleted():
    entry = taxonomy.TaxonomyEntry.build_deleted("2")
    assert entry.tax_id == 2
    assert entry.is_deleted is True
    assert entry.name == "deleted taxid 2"
    assert entry.lineage == ""
    assert entry.replaced_by is None
    row = list(entry.writeable())[0]
    assert row[7] is True


def test_parse_integrates_rank_and_proteome():
    entries = list(
        taxonomy.parse(
            StringIO(LINEAGE_DMP),
            StringIO(NAMES_DMP),
            StringIO(MERGED_DMP),
            StringIO(NODES_DMP),
            StringIO(DELNODES_DMP_EMPTY),
            ref_proteomes_handle=StringIO(STATS_HEADER + STATS_DATA),
        )
    )
    by_id = {e.tax_id: e for e in entries}
    assert by_id[9606].rank == "species"
    assert by_id[9606].reference_proteome is True
    assert by_id[9443].rank == "order"
    assert by_id[9443].reference_proteome is False


def test_parse_without_ref_proteomes():
    entries = list(
        taxonomy.parse(
            StringIO(LINEAGE_DMP),
            StringIO(NAMES_DMP),
            StringIO(MERGED_DMP),
            StringIO(NODES_DMP),
            StringIO(DELNODES_DMP_EMPTY),
        )
    )
    for entry in entries:
        assert entry.reference_proteome is False


def test_build_from_ena(monkeypatch):
    monkeypatch.setattr(
        taxonomy.phy,
        "phylogeny",
        lambda tax_id: {
            "taxId": str(tax_id),
            "scientificName": "uncultured Caudovirales phage",
            "commonName": "",
            "rank": "species",
            "lineage": "Viruses; Duplodnaviria; ",
        },
    )
    entry = taxonomy.TaxonomyEntry.build_from_ena(3016054)
    assert entry.tax_id == 3016054
    assert entry.name == "uncultured Caudovirales phage"
    assert entry.lineage == "Viruses; Duplodnaviria; uncultured Caudovirales phage"
    assert entry.rank == "species"
    assert entry.is_deleted is False
    assert entry.replaced_by is None


def test_resolve_missing_collects_all_failures(monkeypatch):
    def fake_phylogeny(tax_id):
        if tax_id == 9606:
            return {"scientificName": "Homo sapiens", "rank": "species", "lineage": ""}
        if tax_id == 111:
            raise phy.UnknownTaxonId(tax_id)
        raise phy.FailedTaxonId(tax_id)

    monkeypatch.setattr(taxonomy.phy, "phylogeny", fake_phylogeny)

    with pytest.raises(taxonomy.UnresolvableTaxids) as excinfo:
        taxonomy.resolve_missing(["9606", "111", "222"])
    failed_ids = {tid for tid, _ in excinfo.value.failures}
    assert failed_ids == {"111", "222"}


def test_write_missing_round_trips(monkeypatch):
    monkeypatch.setattr(
        taxonomy.phy,
        "phylogeny",
        lambda tax_id: {
            "scientificName": "Homo sapiens",
            "rank": "species",
            "lineage": "Eukaryota; Metazoa; ",
        },
    )
    output = StringIO()
    taxonomy.write_missing(StringIO("9606\n\n"), output)
    row = output.getvalue().strip().split(",")
    assert row[0] == "9606"
    assert row[1] == "Homo sapiens"


def test_parse_emits_deleted_taxids():
    entries = list(
        taxonomy.parse(
            StringIO(LINEAGE_DMP),
            StringIO(NAMES_DMP),
            StringIO(MERGED_DMP),
            StringIO(NODES_DMP),
            StringIO(DELNODES_DMP),
        )
    )
    by_id = {e.tax_id: e for e in entries}
    # 2 is deleted and should be emitted as a deleted sentinel
    assert by_id[2].is_deleted is True
    assert by_id[2].name == "deleted taxid 2"
    # 9606 is in delnodes but still a live node, so it must keep its real data
    assert by_id[9606].is_deleted is False
    assert by_id[9606].name == "Homo sapiens"
    # live taxa are flagged not-deleted
    assert by_id[9443].is_deleted is False
