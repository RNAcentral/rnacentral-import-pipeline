# -*- coding: utf-8 -*-

import io
import json

import pytest

from rnacentral_pipeline.rnacentral.ftp_export import ensembl


def record(**overrides):
    data = {
        "rna_type": "miRNA",
        "sequence": "acgu",
        "xrefs": [
            {
                "database": "miRBase",
                "external_id": "MI1",
                "optional_id": None,
                "molecule_type": None,
            }
        ],
    }
    data.update(overrides)
    return data


def lines(*records):
    return io.StringIO("".join(json.dumps(r) + "\n" for r in records))


def test_modomics_ids_become_urls():
    xref = {"database": "Modomics", "external_id": "12", "optional_id": None}
    assert (
        ensembl.external_id(xref) == "http://modomics.genesilico.pl/sequences/list/12"
    )


def test_disallowed_rna_types_are_not_sent():
    assert ensembl.builder(record(rna_type="circRNA")) == {}


def test_generate_file_writes_the_cleaned_records_as_json():
    out = io.StringIO()
    ensembl.generate_file(lines(record()), out)
    assert json.loads(out.getvalue()) == [
        {
            "rna_type": "miRNA",
            "sequence": "ACGT",
            "xrefs": [{"database": "miRBase", "id": "MI1"}],
        }
    ]


def test_generate_file_drops_disallowed_types_and_invalid_sequences():
    out = io.StringIO()
    ensembl.generate_file(
        lines(record(rna_type="sgRNA"), record(sequence="ACGZ"), record(sequence="")),
        out,
    )
    assert out.getvalue() == ""


def test_generate_file_validates_against_the_schema(tmp_path):
    schema = tmp_path / "schema.json"
    schema.write_text(json.dumps({"type": "array", "items": {"required": ["gene"]}}))
    with pytest.raises(Exception):
        ensembl.generate_file(lines(record()), io.StringIO(), str(schema))
    schema.write_text(
        json.dumps({"type": "array", "items": {"required": ["sequence"]}})
    )
    out = io.StringIO()
    ensembl.generate_file(lines(record()), out, str(schema))
    assert json.loads(out.getvalue())[0]["sequence"] == "ACGT"
