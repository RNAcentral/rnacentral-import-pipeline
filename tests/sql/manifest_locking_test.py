import threading

import psycopg2
import psycopg2.extensions
import pytest

from rnacentral_pipeline.databases import manifest

pytestmark = [pytest.mark.db]


class PausedCommit(psycopg2.extensions.connection):
    # Holds its transaction open at commit until the test lets it go

    paused = threading.Event()
    resume = threading.Event()
    armed = False

    def commit(self):
        if self.armed:
            self.paused.set()
            self.resume.wait(10)
        super().commit()


@pytest.fixture
def pg(pg_url):
    conn = psycopg2.connect(pg_url)
    with conn.cursor() as cur:
        cur.execute("CREATE SCHEMA IF NOT EXISTS rnacen")
    conn.commit()
    manifest.ensure_partition(conn, "ENSEMBL")
    yield pg_url
    with conn.cursor() as cur:
        cur.execute("DROP SCHEMA rnacen CASCADE")
    conn.commit()
    conn.close()


# Both run straight after a release; on release 29 they deadlocked on Ensembl's rows
def test_applying_a_manifest_waits_for_the_sources_of_the_same_database(
    pg, tmp_path, monkeypatch
):
    # The tables exist; its own commit would otherwise be the one paused, before any lock
    monkeypatch.setattr(manifest, "ensure_table", lambda conn: None)
    rows = tmp_path / "manifest.csv"
    rows.write_text("ENSEMBL,ENST1,sig,/ensembl/new.gff3\n")

    sources = psycopg2.connect(pg, connection_factory=PausedCommit)
    sources.armed = True
    writer = threading.Thread(
        target=manifest.store_file_signatures,
        args=(sources, "ENSEMBL", {"/ensembl/old.gff3": "s1"}),
    )
    writer.start()
    assert PausedCommit.paused.wait(10)

    manifests = psycopg2.connect(pg)
    applier = threading.Thread(target=manifest.apply_artifacts, args=(manifests, rows))
    applier.start()
    applier.join(2)
    blocked = applier.is_alive()

    PausedCommit.resume.set()
    writer.join(10)
    applier.join(10)
    sources.close()
    manifests.close()
    assert blocked
