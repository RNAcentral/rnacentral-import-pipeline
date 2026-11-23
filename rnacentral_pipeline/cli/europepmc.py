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

from rnacentral_pipeline.databases.europepmc import stream
from rnacentral_pipeline.output_format import format_option


@click.group("europepmc")
def cli():
    """
    This group of commands deal with handling publication data.
    """
    pass


@cli.command("stream-lookup")
@click.option("--column", default=0)
@click.option("--allow-fallback/--no-allow-fallback", default=False)
@click.option("--ignore-missing/--no-ignore-missing", default=True)
@click.argument("directory", default="out", type=click.Path())
@click.argument("ids", default="ref_ids.csv", type=click.File("r"))
@click.argument("output", default="references.csv", type=click.Path())
@format_option
def stream_lookup(
    directory, ids, output, column=0, allow_fallback=False, ignore_missing=True
):
    """
    Load all ids to write in. Format is governed by
    ``--format``/``RNAC_OUTPUT_FORMAT`` (CSV by default).
    """
    stream.write_lookup(
        ids,
        directory,
        output,
        column=column,
        ignore_missing=ignore_missing,
        allow_fallback=allow_fallback,
    )
