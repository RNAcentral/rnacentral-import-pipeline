"""
Copyright [2009-2020] EMBL-European Bioinformatics Institute
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

from rnacentral_pipeline.databases.rediportal import parser
from rnacentral_pipeline.rnacentral.notify.slack import send_notification


@click.group("rediportal")
def cli():
    """
    A group of commands dealing with RediPortal data.
    """
    pass


@cli.command("parse-bed")
@click.argument("redi_bedfile", type=click.File("r"))
@click.argument("redi_metadata", type=click.File("r"))
@click.argument("rnc_bedfile", type=click.File("r"))
@click.argument("output", type=click.Path())
def parse_rediportal(redi_bedfile, redi_metadata, rnc_bedfile, output):
    """
    Intersect REDIportal bedfile with ours, parse the result alongside the metadata

    """
    parser.parse(redi_bedfile, redi_metadata, rnc_bedfile, output)
