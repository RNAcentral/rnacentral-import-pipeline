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
from . import external
from . import ftp_export
from . import misc
from . import ontologies
from . import pdb
from . import precompute
from . import rfam


@click.group()
def cli():
    """
    This script contains commands for dealing with the RNAcentral import
    pipeline. This handles individual python parts, and the overall pipeline is
    coordinated by nextflow.
    """
    pass


cli.add_command(ensembl.cli)
cli.add_command(external.cli)
cli.add_command(ftp_export.cli)
cli.add_command(ontologies.cli)
cli.add_command(precompute.cli)
cli.add_command(rfam.cli)
cli.add_command(misc.run_release)
cli.add_command(misc.find_upi_ranges)
cli.add_command(misc.crs_data)
cli.add_command(pdb.cli)