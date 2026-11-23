# -*- coding: utf-8 -*-

import io
import json

from rnacentral_pipeline.rnacentral.search_export import metadata


def lines(*rows):
    return io.StringIO("".join(json.dumps(r) + "\n" for r in rows))


def test_so_term_tree_writes_one_json_line_per_type(monkeypatch):
    monkeypatch.setattr(metadata.so, "load_ontology", lambda url: "ontology")
    monkeypatch.setattr(metadata.so, "rna_type_tree", lambda ont, term: [[term, ont]])
    out = io.StringIO()

    metadata.write_so_term_tree(io.StringIO("SO:1\nSO:2\n"), "url", out)

    assert [json.loads(line) for line in out.getvalue().splitlines()] == [
        {"so_rna_type": "SO:1", "so_term_tree": [["SO:1", "ontology"]]},
        {"so_rna_type": "SO:2", "so_term_tree": [["SO:2", "ontology"]]},
    ]
