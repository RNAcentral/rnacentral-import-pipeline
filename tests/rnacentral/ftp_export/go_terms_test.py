# -*- coding: utf-8 -*-

import io
import json

from rnacentral_pipeline.rnacentral.ftp_export import go_terms


def line(taxid, **extra):
    return json.dumps(
        {
            "id": "RF1",
            "ontology_term_id": "GO:1",
            "models": "m",
            "taxid": taxid,
            **extra,
        }
    )


def test_export_drops_uncultured_organisms():
    handle = io.StringIO("\n".join([line(9606), line(155900)]))
    out = io.StringIO()
    go_terms.export(handle, out)
    assert out.getvalue().splitlines() == ["RF1\tGO:1\tm"]


def test_export_ignores_extra_fields_and_accepts_string_taxids():
    out = io.StringIO()
    go_terms.export(io.StringIO(line("198431") + "\n" + line("1", other="x")), out)
    assert out.getvalue().splitlines() == ["RF1\tGO:1\tm"]
