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

from rnacentral_pipeline.databases.mgnify import prepare

# ``inferredPhylogeny`` with none of the d__/p__/.../s__ markers, so every
# gtdb helper returns None immediately without needing a taxid lookup
# (network/cassette), and prepare_mgnify_data falls back to its own
# catalog_name -> taxid table.
NO_PHYLOGENY = "unclassified"


def raw_data(*entries):
    return {"metaData": {"source": "test"}, "data": list(entries)}


def entry(source_model, catalog_name="human gut genome catalogue"):
    return {
        "inferredPhylogeny": NO_PHYLOGENY,
        "additionalAnnotations": {"catalog_name": catalog_name},
        "sourceModel": source_model,
    }


def test_uses_trnascanse_lookup_directly(monkeypatch):
    monkeypatch.setattr(prepare, "get_so_type", lambda conn_str: {})

    result = prepare.prepare_mgnify_data(raw_data(entry("TRNASCANSE:Ala")), "conn_str")

    (prepared,) = result["data"]
    assert prepared["soTermId"] == "SO:0000254"


def test_uses_rfam_model_lookup_for_other_sources(monkeypatch):
    monkeypatch.setattr(
        prepare, "get_so_type", lambda conn_str: {"my_model": "SO:0000655"}
    )

    result = prepare.prepare_mgnify_data(raw_data(entry("RF:my_model")), "conn_str")

    (prepared,) = result["data"]
    assert prepared["soTermId"] == "SO:0000655"


def test_falls_back_to_catalog_taxid_when_phylogeny_is_unresolved(monkeypatch):
    monkeypatch.setattr(prepare, "get_so_type", lambda conn_str: {})

    result = prepare.prepare_mgnify_data(
        raw_data(entry("TRNASCANSE:Ala", catalog_name="human gut genome catalogue")),
        "conn_str",
    )

    (prepared,) = result["data"]
    assert prepared["taxonId"] == "NCBITaxon:408170"


def test_preserves_metadata_and_processes_every_entry(monkeypatch):
    monkeypatch.setattr(prepare, "get_so_type", lambda conn_str: {})

    result = prepare.prepare_mgnify_data(
        raw_data(entry("TRNASCANSE:Ala"), entry("TRNASCANSE:Val")), "conn_str"
    )

    assert result["metaData"] == {"source": "test"}
    assert len(result["data"]) == 2
