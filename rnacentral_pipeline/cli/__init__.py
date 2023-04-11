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

from rnacentral_pipeline.cli import (
    context,
    cpat,
    crw,
    ena,
    ensembl,
    europepmc,
    evlncrnas,
    expressionatlas,
    five_s_rrnadb,
    flybase,
    ftp_export,
    genecards_suite,
    genes,
    genome_mapping,
    gtrnadb,
    hgnc,
    intact,
    lncbase,
    lncbook,
    lncipedia,
    mirbase,
    mirgenedb,
    misc,
    ncbi,
    notify,
    ols,
    pdb,
    pirbase,
    plncdb,
    pombase,
    precompute,
    psicquic,
    qa,
    quickgo,
    r2dt,
    refseq,
    release,
    repeats,
    rfam,
    ribocentre,
    ribovision,
    scan_imports,
    search_export,
    sgd,
    silva,
    snodb,
    snorna_database,
    tarbase,
    zfin,
    zwd,
)


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


cli.add_command(context.cli)
cli.add_command(cpat.cli)
cli.add_command(crw.cli)
cli.add_command(ena.cli)
cli.add_command(ensembl.cli)
cli.add_command(europepmc.cli)
cli.add_command(evlncrnas.cli)
cli.add_command(expressionatlas.cli)
cli.add_command(five_s_rrnadb.cli)
cli.add_command(flybase.cli)
cli.add_command(ftp_export.cli)
cli.add_command(genecards_suite.cli)
cli.add_command(genes.cli)
cli.add_command(genome_mapping.cli)
cli.add_command(gtrnadb.cli)
cli.add_command(hgnc.cli)
cli.add_command(intact.cli)
cli.add_command(lncbase.cli)
cli.add_command(lncbook.cli)
cli.add_command(lncipedia.cli)
cli.add_command(mirbase.cli)
cli.add_command(mirgenedb.cli)
cli.add_command(misc.crs_data)
cli.add_command(misc.find_upi_ranges)
cli.add_command(misc.validate_pgloader)
cli.add_command(ncbi.cli)
cli.add_command(notify.cli)
cli.add_command(ols.cli)
cli.add_command(pdb.cli)
cli.add_command(pirbase.cli)
cli.add_command(plncdb.cli)
cli.add_command(pombase.cli)
cli.add_command(psicquic.cli)
cli.add_command(precompute.cli)
cli.add_command(qa.cli)
cli.add_command(quickgo.cli)
cli.add_command(r2dt.cli)
cli.add_command(refseq.cli)
cli.add_command(release.cli)
cli.add_command(repeats.cli)
cli.add_command(rfam.cli)
cli.add_command(ribovision.cli)
cli.add_command(ribocentre.cli)
cli.add_command(scan_imports.cli)
cli.add_command(search_export.cli)
cli.add_command(sgd.cli)
cli.add_command(silva.cli)
cli.add_command(snodb.cli)
cli.add_command(snorna_database.cli)
cli.add_command(tarbase.cli)
cli.add_command(zfin.cli)
cli.add_command(zwd.cli)
