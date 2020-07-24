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

import click

from . import ensembl
from . import europepmc
from . import external
from . import fetch
from . import ftp_export
from . import genome_mapping
from . import lookup
from . import misc
from . import ncbi
from . import ols
from . import precompute
from . import qa
from . import rfam
from . import search_export
from . import text_mining
from . import r2dt


@click.group()
def cli():
    """
    This script contains commands for dealing with the RNAcentral import
    pipeline. This handles individual python parts, and the overall pipeline is
    coordinated by nextflow.
    """
    pass


cli.add_command(ensembl.cli)
cli.add_command(europepmc.cli)
cli.add_command(external.cli)
cli.add_command(fetch.cli)
cli.add_command(ftp_export.cli)
cli.add_command(genome_mapping.cli)
cli.add_command(lookup.cli)
cli.add_command(misc.check_release)
cli.add_command(misc.crs_data)
cli.add_command(misc.find_upi_ranges)
cli.add_command(misc.run_release)
cli.add_command(misc.validate_pgloader)
cli.add_command(ncbi.cli)
cli.add_command(ols.cli)
cli.add_command(precompute.cli)
cli.add_command(qa.cli)
cli.add_command(rfam.cli)
cli.add_command(search_export.cli)
cli.add_command(text_mining.cli)
cli.add_command(r2dt.cli)
