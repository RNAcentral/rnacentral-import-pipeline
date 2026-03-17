# -*- coding: utf-8 -*-

import json
from pathlib import Path

import pandas as pd

from rnacentral_pipeline.databases.evlncrnas import helpers, parser


def test_can_parse_evlncrnas3_fixture(monkeypatch, tmp_path):
    fixture_dir = Path("data/evlncrnas")

    monkeypatch.setattr(
        parser,
        "handled_phylogeny",
        lambda species: {
            "Homo sapiens": 9606,
        }.get(species),
    )

    def enrich(frame):
        enriched = frame.copy()
        enriched["sequence"] = "ACGTACGTACGT"
        enriched["assembly_id"] = "GRCh38"
        enriched["chromosome"] = "11"
        enriched["region_start"] = 65497688
        enriched["region_stop"] = 65506516
        enriched["chain"] = "+"
        return enriched, enriched.iloc[0:0].copy()

    def fake_empty_matches(_frame, _dump):
        return pd.DataFrame(
            columns=[
                "ID",
                "external_id",
                "lookup_name",
                "urs",
                "taxid",
                "is_exact_match",
            ]
        )

    monkeypatch.setattr(parser, "get_ensembl_accessions", enrich)
    monkeypatch.setattr(parser, "get_ncbi_accessions", enrich)
    monkeypatch.setattr(parser, "get_db_matches", fake_empty_matches)
    monkeypatch.setattr(parser.lookup, "as_mapping", lambda *_args, **_kwargs: {})
    monkeypatch.setattr(helpers, "lineage", lambda _record: "Eukaryota; Metazoa; Mammalia")

    dump = tmp_path / "ev_lookup.csv"
    dump.write_text("urs,taxid,external_id\n", encoding="utf-8")

    entries = list(parser.parse(fixture_dir, (dump,), "postgres://ignored"))

    assert len(entries) == 1

    entry = entries[0]
    assert entry.primary_id == "EVLNCRNAS:EL3692"
    assert entry.accession == "EVLNCRNAS:EL3692"
    assert entry.ncbi_tax_id == 9606
    assert entry.database == "EVLNCRNAS"
    assert entry.sequence == "ACGTACGTACGT"
    assert entry.rna_type == "SO:0001463"
    assert entry.url == "https://www.sdklab-biophysics-dzu.net/EVLncRNAs3/#/detail?id=EL3692"
    assert entry.gene_synonyms == [
        "MALAT1",
        "HCN",
        "LINC00047",
        "NCRNA00047",
        "NEAT2",
        "PRO2853",
    ]
    assert sorted(ref.external_id for ref in entry.references) == ["31174563", "33235630"]
    assert entry.note_data == {}
    assert json.loads(entry.note) == {
        "url": "https://www.sdklab-biophysics-dzu.net/EVLncRNAs3/#/detail?id=EL3692"
    }
