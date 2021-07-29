# -*- coding: utf-8 -*-

"""
Copyright [2009-2021] EMBL-European Bioinformatics Institute
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

from pathlib import Path
import pickle

import click

from rnacentral_pipeline.cpat import parser
from rnacentral_pipeline.cpat.data import CpatWriter
from rnacentral_pipeline.writers import build


@click.group("cpat")
def cli():
    """
    Commands with processing the Rfam metadata.
    """
    pass


@cli.command("parse")
@click.argument("cutoffs", type=click.File("rb"))
@click.argument("model_name")
@click.argument("results", type=click.Path())
@click.argument("output", type=click.Path())
def parse(cutoffs, model_name, results, output):
    cutoffs = pickle.load(cutoffs)
    data = parser.parse(cutoffs, model_name, Path(results))
    with build(CpatWriter, Path(output)) as wtr:
        wtr.write(data)


@cli.command("generate-cutoffs")
@click.argument("data-folder", type=click.Path())
@click.argument("output", type=click.File("wb"))
def generate_cutoffs(data_folder, output):
    cutoffs = parser.cutoffs(Path(data_folder))
    pickle.dump(cutoffs, output)
