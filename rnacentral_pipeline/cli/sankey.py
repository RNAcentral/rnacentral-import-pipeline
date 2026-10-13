# -*- coding: utf-8 -*-

from pathlib import Path

import click

from rnacentral_pipeline.rnacentral import sankey


@click.command("sankey")
@click.argument("output", type=click.Path(path_type=Path))
@click.argument("names", nargs=-1)
@click.option("--db-url", envvar="PGDATABASE")
def cli(output: Path, names, db_url: str):
    """
    Draw a sankey plot into OUTPUT for each database in NAMES, or each database
    loaded since its plot was last drawn, plus one for all of them together.
    """
    output.mkdir(parents=True, exist_ok=True)
    sankey.write_all(db_url, output, names)
