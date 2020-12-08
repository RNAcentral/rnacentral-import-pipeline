# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import logging

import click

from . import (ena, ensembl, europepmc, external, ftp_export, genes,
               genome_mapping, gtrnadb, lookup, misc, ncbi, ols, pdb, pirbase,
               precompute, qa, r2dt, rfam, search_export, text_mining, zfin)


@click.group()
@click.option(
    "--log-level",
    default="warn",
    type=click.Choice(
        ["critical", "error", "warn", "info", "debug"], case_sensitive=False
    ),
)
def cli(log_level):
    """
    This script contains commands for dealing with the RNAcentral import
    pipeline. This handles individual python parts, and the overall pipeline is
    coordinated by nextflow.
    """
    level = getattr(logging, log_level.upper())
    logger = logging.getLogger()
    logger.setLevel(level=level)
    pass


cli.add_command(ena.cli)
cli.add_command(ensembl.cli)
cli.add_command(europepmc.cli)
cli.add_command(external.cli)
cli.add_command(ftp_export.cli)
cli.add_command(genes.cli)
cli.add_command(genome_mapping.cli)
cli.add_command(gtrnadb.cli)
cli.add_command(lookup.cli)
cli.add_command(misc.check_release)
cli.add_command(misc.crs_data)
cli.add_command(misc.find_upi_ranges)
cli.add_command(misc.run_release)
cli.add_command(misc.validate_pgloader)
cli.add_command(ncbi.cli)
cli.add_command(ols.cli)
cli.add_command(pirbase.cli)
cli.add_command(pdb.cli)
cli.add_command(precompute.cli)
cli.add_command(qa.cli)
cli.add_command(r2dt.cli)
cli.add_command(rfam.cli)
cli.add_command(search_export.cli)
cli.add_command(text_mining.cli)
cli.add_command(zfin.cli)
