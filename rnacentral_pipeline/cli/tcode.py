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

import os
from pathlib import Path

import click

from rnacentral_pipeline import schemas
from rnacentral_pipeline.parquet_writers import typed_parquet_writer
from rnacentral_pipeline.tcode import parser
from rnacentral_pipeline.tcode.data import TcodeWriter
from rnacentral_pipeline.writers import build


@click.group("tcode")
def cli():
    """
    Commands for parsing TCODE output.
    """
    pass


@cli.command("parse")
@click.argument("tcode_output", type=click.Path())
@click.argument("output", type=click.Path())
def parse(tcode_output, output):
    data = parser.parse(Path(tcode_output))
    out_path = Path(output)
    if os.environ.get("RNAC_OUTPUT_FORMAT", "csv").lower() == "parquet":
        out_path.mkdir(parents=True, exist_ok=True)
        with typed_parquet_writer(
            out_path / "results.parquet", schemas.TCODE_RESULTS
        ) as writer:
            for result in data:
                writer.writerow(result.writeable())
        return
    with build(TcodeWriter, out_path) as wtr:
        wtr.write(data)
