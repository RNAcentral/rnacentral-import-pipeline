# -*- coding: utf-8 -*-

from pathlib import Path

import click

from rnacentral_pipeline.stopfree import scan


@click.group("stopfree")
def cli():
    """
    Commands for running stopfree.
    """
    pass


@cli.command("scan")
@click.argument("input_fasta", type=click.Path(exists=True, path_type=Path))
@click.argument("output", type=click.Path(path_type=Path))
@click.option(
    "--max-probability",
    default=0.05,
    show_default=True,
    type=float,
    help="Maximum null-model probability still considered protein coding.",
)
def scan_command(input_fasta: Path, output: Path, max_probability: float):
    scan.run(input_fasta, output, max_probability)
