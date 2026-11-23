# -*- coding: utf-8 -*-

import os
import shutil
import subprocess
from pathlib import Path

import pytest

RELATED = Path("files/import-data/post-release/001__related-sequences.sql")

pytestmark = [pytest.mark.db]


def run_related(pg_url, path=RELATED):
    psql = shutil.which("psql")
    if not psql:
        pytest.skip("No psql")
    env = dict(os.environ, PGOPTIONS="-c search_path=public,rnacen")
    return subprocess.run(
        [psql, "-X", "-q", "-v", "ON_ERROR_STOP=1", "-f", str(path), pg_url],
        capture_output=True,
        text=True,
        env=env,
    )


def seed(db):
    with db.cursor() as cur:
        cur.execute("INSERT INTO rnc_database VALUES (25, 'ENSEMBL')")
        cur.execute(
            "INSERT INTO rnc_accessions VALUES ('ENST_S', 'ENSEMBL', NULL), "
            "('ENST_T', 'ENSEMBL', NULL), ('ENST_G', 'ENSEMBL', 'ENSG1.3'), "
            "('ENST_OLD', 'ENSEMBL', NULL)"
        )
        cur.execute(
            "INSERT INTO xref (ac, urs, taxid, deleted, dbid) VALUES "
            "('ENST_S', 'URS1', 9606, 'N', 25), ('ENST_T', 'URS2', 9606, 'N', 25), "
            "('ENST_G', 'URS3', 9606, 'N', 25), ('ENST_OLD', 'URS5', 9606, 'Y', 25)"
        )
        cur.execute(
            "INSERT INTO rnacen.rnc_related_sequences VALUES "
            "('URS1_9606', 'ENST_S', 'URS2_9606', 'ENST_T', 'isoform', '{old}')"
        )
        cur.executemany(
            "INSERT INTO load_rnc_related_sequences "
            "(source_accession, target_accession, relationship_type, methods) "
            "VALUES ('ENST_S', %s, %s, '{m1}')",
            [
                ("ENST_T", "isoform"),
                ("ENST_T", "isoform"),
                ("ENST_OLD", "isoform"),
                ("UNKNOWN", "target_rna"),
                ("ENSEMBL:ENSG1", "target_rna"),
                ("P12345", "target_protein"),
            ],
        )


def related(db):
    with db.cursor() as cur:
        cur.execute(
            "SELECT target_accession, relationship_type::text, source_urs_taxid, "
            "target_urs_taxid, ARRAY(SELECT unnest(methods) ORDER BY 1) "
            "FROM rnacen.rnc_related_sequences ORDER BY 1, 2"
        )
        return cur.fetchall()


def test_targets_resolve_to_their_urs(db, pg_url):
    seed(db)
    result = run_related(pg_url)
    assert result.returncode == 0, result.stderr
    assert related(db) == [
        ("ENSEMBL:ENSG1", "target_rna", "URS1_9606", "URS3_9606", ["m1"]),
        ("ENST_OLD", "isoform", "URS1_9606", None, ["m1"]),
        ("ENST_T", "isoform", "URS1_9606", "URS2_9606", ["m1", "old"]),
        ("P12345", "target_protein", "URS1_9606", None, ["m1"]),
        ("UNKNOWN", "target_rna", "URS1_9606", None, ["m1"]),
    ]
