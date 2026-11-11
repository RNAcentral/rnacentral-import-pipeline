# -*- coding: utf-8 -*-

import csv

from click.testing import CliRunner

from rnacentral_pipeline.cli import ensembl as ensembl_cli
from rnacentral_pipeline.databases import manifest
from rnacentral_pipeline.databases.ensembl import delta

EMBL = "http://ftp/Homo_sapiens/ensembl/geneset/2026_04/genes.embl.gz"
GFF = "http://ftp/Homo_sapiens/ensembl/geneset/2026_04/genes.gff3.gz"
MD5S = "aaa  ./genes.embl.gz\nbbb  ./genes.gff3.gz\nccc  ./genes.gtf.gz\n"


def row(name, url):
    return (name, "9606", url, url.replace("embl", "gff3"), "GCA_1")


def test_signature_follows_the_checksums_of_the_two_files_we_read():
    assert delta.source_signature(
        EMBL, GFF, lambda url: MD5S
    ) == delta.source_signature(EMBL, GFF, lambda url: MD5S)
    other = delta.source_signature(EMBL, GFF, lambda url: MD5S.replace("bbb", "zzz"))
    assert other != delta.source_signature(EMBL, GFF, lambda url: MD5S)


def test_signature_ignores_files_we_do_not_read():
    other = MD5S.replace("ccc", "zzz")
    assert delta.source_signature(
        EMBL, GFF, lambda url: MD5S
    ) == delta.source_signature(EMBL, GFF, lambda url: other)


def test_signature_reads_the_listing_from_the_geneset_directory():
    seen = []
    delta.source_signature(EMBL, GFF, lambda url: seen.append(url) or MD5S)
    assert seen == ["http://ftp/Homo_sapiens/ensembl/geneset/2026_04/md5sum.txt"]


def test_signature_never_matches_when_the_checksums_cannot_be_read(monkeypatch):
    def refuse(url):
        raise OSError("refused")

    monkeypatch.setattr(delta.time, "time", iter([1.0, 2.0]).__next__)
    assert delta.source_signature(EMBL, GFF, refuse) != delta.source_signature(
        EMBL, GFF, refuse
    )


def test_signature_never_matches_when_a_file_is_missing_from_the_listing():
    assert delta.source_signature(
        EMBL, GFF, lambda url: "aaa  ./genes.embl.gz\n"
    ).startswith("unverified")


def test_signatures_read_the_listing_from_the_nfs_copy(tmp_path):
    geneset = tmp_path / "Homo_sapiens/ensembl/geneset/2026_04"
    geneset.mkdir(parents=True)
    (geneset / "md5sum.txt").write_text(MD5S)
    read = delta.mirror("http://ftp", str(tmp_path))
    assert delta.signatures([("h", "9606", EMBL, GFF, "GCA_1")], read) == {
        EMBL: delta.source_signature(EMBL, GFF, lambda url: MD5S)
    }


def test_only_new_and_changed_genesets_are_fetched():
    rows = [
        row("a", "u/a/genes.embl.gz"),
        row("b", "u/b/genes.embl.gz"),
        row("c", "u/c/genes.embl.gz"),
    ]
    current = {
        "u/a/genes.embl.gz": "1",
        "u/b/genes.embl.gz": "2",
        "u/c/genes.embl.gz": "3",
    }
    stored = {"u/a/genes.embl.gz": "1", "u/b/genes.embl.gz": "old"}
    fetch, gone, sources = delta.plan(rows, rows, stored, current)
    assert [r[0] for r in fetch] == ["b", "c"]
    assert fetch[0][-1] == "2"
    assert gone == []
    assert sources == current


def test_a_geneset_that_left_the_list_is_gone():
    rows = [row("a", "u/a/genes.embl.gz")]
    stored = {"u/a/genes.embl.gz": "1", "u/old/genes.embl.gz": "9"}
    fetch, gone, sources = delta.plan(rows, rows, stored, {"u/a/genes.embl.gz": "1"})
    assert gone == ["u/old/genes.embl.gz"]
    assert "u/old/genes.embl.gz" not in sources


# A run limited to some species must neither retire the rest nor mark them up to date.
def test_species_left_out_of_a_partial_run_are_not_gone_and_keep_their_signature():
    rows = [
        row("a", "u/a/genes.embl.gz"),
        row("b", "u/b/genes.embl.gz"),
        row("c", "u/c/genes.embl.gz"),
    ]
    picked = delta.select(rows, ["A"])
    assert [r[0] for r in picked] == ["a"]
    stored = {"u/b/genes.embl.gz": "old-b"}
    fetch, gone, sources = delta.plan(rows, picked, stored, {"u/a/genes.embl.gz": "1"})
    assert [r[0] for r in fetch] == ["a"]
    assert gone == []
    assert sources == {"u/a/genes.embl.gz": "1", "u/b/genes.embl.gz": "old-b"}


def write_manifest(path, rows):
    with path.open("w", newline="") as handle:
        csv.writer(handle).writerows(rows)


def stored_dump(rows):
    def dump(conn, database, path):
        with path.open("w", newline="") as handle:
            csv.writer(handle).writerows(rows.get(database, []))

    return dump


def test_accessions_a_reparsed_geneset_lost_are_deleted(tmp_path, monkeypatch):
    monkeypatch.setattr(
        manifest,
        "dump_signatures",
        stored_dump(
            {"ENSEMBL": [("T1", "s", "u/a"), ("T2", "s", "u/a"), ("T3", "s", "u/b")]}
        ),
    )
    write_manifest(tmp_path / "m.csv", [("ENSEMBL", "T1", "s2", "u/a")])
    found = delta.deletions(None, tmp_path / "m.csv", ["u/a"], tmp_path)
    assert found == [("ENSEMBL", "T2")]


def test_a_gone_geneset_loses_all_its_accessions(tmp_path, monkeypatch):
    monkeypatch.setattr(
        manifest, "dump_signatures", stored_dump({"TAIR": [("TAIR:AT1", "s", "u/a")]})
    )
    write_manifest(tmp_path / "m.csv", [])
    found = delta.deletions(None, tmp_path / "m.csv", ["u/a"], tmp_path)
    assert found == [("TAIR", "TAIR:AT1")]


def test_the_first_run_deletes_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(manifest, "dump_signatures", stored_dump({}))
    write_manifest(tmp_path / "m.csv", [("ENSEMBL", "T1", "s", "u/a")])
    assert delta.deletions(None, tmp_path / "m.csv", ["u/a"], tmp_path) == []


def run_source_diff(tmp_path, monkeypatch, stored, *options):
    listing = tmp_path / "species.csv"
    write_manifest(
        listing, [row("a", "u/a/genes.embl.gz"), row("b", "u/b/genes.embl.gz")]
    )
    monkeypatch.setattr(
        ensembl_cli.delta, "signatures", lambda rows, read: {r[2]: "new" for r in rows}
    )
    monkeypatch.setattr(
        ensembl_cli.manifest, "load_file_signatures_for", lambda url, db: dict(stored)
    )
    names = ("to_fetch.csv", "scanned.txt", "sources.csv")
    result = CliRunner().invoke(
        ensembl_cli.cli,
        ["source-diff", "--ftp-url", "u", "--nfs-path", "n", *options, str(listing)]
        + [str(tmp_path / n) for n in names],
    )
    assert result.exit_code == 0, result.output
    return [(tmp_path / n).read_text() for n in names]


def test_source_diff_scans_the_changed_and_the_gone(tmp_path, monkeypatch):
    stored = {
        "u/a/genes.embl.gz": "new",
        "u/b/genes.embl.gz": "old",
        "u/x/genes.embl.gz": "x",
    }
    to_fetch, scanned, sources = run_source_diff(tmp_path, monkeypatch, stored)
    assert [r[0] for r in csv.reader(to_fetch.splitlines())] == ["b"]
    assert scanned.split() == ["u/b/genes.embl.gz", "u/x/genes.embl.gz"]
    assert sources.splitlines() == [
        "ENSEMBL,u/a/genes.embl.gz,new",
        "ENSEMBL,u/b/genes.embl.gz,new",
    ]


def test_force_full_fetches_every_geneset_and_retires_none(tmp_path, monkeypatch):
    stored = {"u/a/genes.embl.gz": "new", "u/x/genes.embl.gz": "x"}
    to_fetch, scanned, _ = run_source_diff(
        tmp_path, monkeypatch, stored, "--force-full"
    )
    assert len(to_fetch.splitlines()) == 2
    assert "u/x/genes.embl.gz" not in scanned
