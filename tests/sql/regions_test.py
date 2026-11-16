# -*- coding: utf-8 -*-

import shutil
import subprocess
from pathlib import Path

import pytest

REGIONS = Path("files/import-data/post-release/001__regions.sql")

pytestmark = [pytest.mark.db]


def run_regions(pg_url):
    psql = shutil.which("psql")
    if not psql:
        pytest.skip("No psql")
    return subprocess.run(
        [psql, "-X", "-q", "-v", "ON_ERROR_STOP=1", "-f", str(REGIONS), pg_url],
        capture_output=True,
        text=True,
    )


def seed(db, load_rows):
    with db.cursor() as cur:
        cur.execute("INSERT INTO ensembl_assembly VALUES ('GRCh38', 'hg38')")
        cur.execute(
            "INSERT INTO xref VALUES ('ENST1', 'URS1', 9606, 'N'), "
            "('ENST2', 'URS2', 9606, 'N'), ('ENST3', 'URS3', 9606, 'N')"
        )
        cur.execute(
            "INSERT INTO rnc_rna_precomputed VALUES ('URS1_9606', false), "
            "('URS2_9606', false), ('URS3_9606', false)"
        )
        for accession, name, identity in [
            ("ENST1", "URS1_9606@1/100-200:+", None),
            ("ENST3", "URS3_9606@1/500-600:+", 0.95),
        ]:
            cur.execute(
                "INSERT INTO rnc_sequence_regions (urs_taxid, region_name, chromosome, "
                "strand, region_start, region_stop, assembly_id, exon_count, was_mapped, identity) "
                "VALUES (%s, %s, '1', 1, 100, 200, 'GRCh38', 1, false, %s) RETURNING id",
                (name.split("@")[0], name, identity),
            )
            region_id = cur.fetchone()[0]
            start = int(name.split("/")[1].split("-")[0])
            cur.execute(
                "INSERT INTO rnc_sequence_exons (region_id, exon_start, exon_stop) "
                "VALUES (%s, %s, %s)",
                (region_id, start, start + 100),
            )
            cur.execute(
                "INSERT INTO rnc_accession_sequence_region VALUES (%s, %s)",
                (accession, region_id),
            )
        cur.executemany(
            "INSERT INTO load_rnc_sequence_regions (accession, region_name, chromosome, "
            "strand, exon_start, exon_stop, assembly_id, exon_count) "
            "VALUES (%s, %s, '1', 1, %s, %s, 'GRCh38', %s)",
            load_rows,
        )


def region(db, name):
    with db.cursor() as cur:
        cur.execute(
            "SELECT xmin::text, identity, "
            "(SELECT count(*) FROM rnc_sequence_exons e WHERE e.region_id = r.id) "
            "FROM rnc_sequence_regions r WHERE region_name = %s",
            (name,),
        )
        return cur.fetchone()


def test_reload_leaves_unchanged_regions_alone(db, pg_url):
    seed(
        db,
        [
            ("ENST1", "@1/100-200:+", 100, 200, 1),
            ("ENST2", "@1/300-400:+", 300, 340, 2),
            ("ENST2", "@1/300-400:+", 360, 400, 2),
            ("ENST3", "@1/500-600:+", 500, 600, 1),
        ],
    )
    unchanged_xmin = region(db, "URS1_9606@1/100-200:+")[0]

    result = run_regions(pg_url)
    assert result.returncode == 0, result.stderr

    assert region(db, "URS1_9606@1/100-200:+")[0] == unchanged_xmin
    assert region(db, "URS2_9606@1/300-400:+")[2] == 2
    assert region(db, "URS3_9606@1/500-600:+")[1] is None
    with db.cursor() as cur:
        cur.execute("SELECT bool_and(has_coordinates) FROM rnc_rna_precomputed")
        assert cur.fetchone()[0] is True


def test_region_missing_exons_fails(db, pg_url):
    seed(
        db,
        [
            ("ENST2", "@1/300-400:+", 300, 340, 3),
            ("ENST2", "@1/300-400:+", 360, 400, 3),
        ],
    )

    result = run_regions(pg_url)

    assert result.returncode != 0
    assert "Some regions 1 are missing exons" in result.stderr
