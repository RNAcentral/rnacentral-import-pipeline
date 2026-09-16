import json

from click.testing import CliRunner

from rnacentral_pipeline.cli import hgnc as hgnc_cli
from rnacentral_pipeline.databases import manifest


def _boom(*args, **kwargs):
    raise AssertionError("HGNC must never parse against the stored manifest")


def test_hgnc_map_never_reads_the_stored_manifest(tmp_path, monkeypatch):
    """
    get_load_release_type pins HGNC to FULL, which retires every xref absent
    from the load. A parse against the manifest would emit only the changed
    records and the release would retire the rest, so the parser must ignore
    the manifest entirely while that pin stands.
    """
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps({"response": {"docs": []}}))
    monkeypatch.setattr(manifest, "load_signatures_for", _boom)
    seen = {}

    def fake_parse(path, db_url, previous):
        seen["previous"] = previous
        return hgnc_cli.parser.ParseResult(
            entries=iter([]), signatures={}, deletions=[]
        )

    monkeypatch.setattr(hgnc_cli.parser, "parse", fake_parse)
    CliRunner().invoke(hgnc_cli.cli, ["map", str(raw), str(tmp_path)])

    assert seen["previous"] == {}
