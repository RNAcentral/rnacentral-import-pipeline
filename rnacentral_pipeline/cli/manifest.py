# -*- coding: utf-8 -*-

"""
CLI for the incremental-parsing manifest (docs/incremental-parsing.md).
"""

import click
import psycopg2

from rnacentral_pipeline.databases import manifest


@click.group("manifest")
def cli():
    """
    Commands for the incremental-parsing import manifest.
    """


@cli.command("apply")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("manifest_csv", type=click.Path(exists=True))
@click.argument("deletions_csv", type=click.Path(), required=False)
def apply(manifest_csv, deletions_csv=None, db_url=None):
    """
    Promote a delta parse's manifest into rnc_import_manifest.

    Run this only AFTER the database's load/release has committed, so a failed load
    leaves the previous manifest untouched and the next run re-emits the same delta.
    """
    conn = psycopg2.connect(db_url)
    try:
        manifest.apply_artifacts(conn, manifest_csv, deletions_csv)
    finally:
        conn.close()
