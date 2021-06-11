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

from pathlib import Path

import click

from rnacentral_pipeline.rnacentral import attempted, r2dt


@click.group("r2dt")
def cli():
    """
    A group of commands for parsing data from secondary structures into an
    importable format.
    """
    pass


@cli.command("process-svgs")
@click.option("--allow-missing", is_flag=True, default=False)
@click.argument("model_info", type=click.File("r"))
@click.argument("directory", type=click.Path())
@click.argument("output", type=click.File("w"))
def process_svgs(model_info, directory, output, allow_missing=False):
    """
    Process all SVG secondary structures in the given directory and produce a
    single data file that can be imported into the database.
    """
    r2dt.write(model_info, directory, output, allow_missing=allow_missing)


@cli.group("should-show")
def should_show():
    """
    Some commands relating to building a model for should show as well as
    running it.
    """


@should_show.command("convert-sheet")
@click.argument("filename", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def convert_sheet(filename, output):
    """
    This command is to convert a downloaded google sheet csv into a csv that can
    be used for training data. Often we will build a spreadsheet of example URS
    and then use that to build a training set. It is nice since you can embedd
    an SVG in google sheets so it is fast for us to compare several of them.

    In order to move that back into the training data you can download that
    sheet as a CSV and then run this command on it to build the CSV that is used
    in training. It requires there be a 'urs' and 'Labeled Should show' column
    to build the CSV. The values in labeled should show must be true/false
    (ignoring case).
    """
    r2dt.write_converted_sheet(filename, output)


@should_show.command("fetch-data")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("filename", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def fetch_training_data(filename, output, db_url=None):
    """
    This builds a CSV file of training data to use for the model building. I
    keep it separate so I can build a training csv and play with it interactivly
    before committing the final modeling building logic to the pipeline.
    """
    r2dt.write_training_data(filename, db_url, output)


@should_show.command("inspect-data")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("filename", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def fetch_inspect_data(filename, output, db_url=None):
    """
    This is the command to use when trying to fetch more examples to add to the
    training set. This will fetch some information that is useful for a person
    to evaluate a diagram and decide if it should be true/false in the training
    set.
    """
    r2dt.write_training_data(filename, db_url, output)


@should_show.command("build-model")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("training-info", type=click.File("r"))
@click.argument("model", type=click.Path())
def build_model(training_info, model, db_url=None):
    """
    This builds a model given then training information. The training
    information should be a csv file of:
        URS,flag
    The flag must be 1 or 0 to indicate if the URS should be shown or not. THis
    will fetch the data like the fetch-data command but will then build a model
    and write it out the the output file directly.
    """
    r2dt.build_model(training_info, db_url, Path(model))


@should_show.command("compute")
@click.option("--db-url", envvar="PGDATABASE")
@click.argument("model", type=click.Path())
@click.argument("filename", type=click.File("r"))
@click.argument("output", type=click.File("w"))
def write_should_show(model, filename, output, db_url=None):
    """
    This computes the should show values for the data in the given file and a
    file listing urs ids to use. The data needed for the URS will be fetched
    from the database. This is meant to operate on large batches, like
    relabeling the entire database.
    """
    r2dt.write_should_show(model, filename, db_url, output)


@cli.group("model-info")
def model_info():
    """
    Commands for parsing and generating data files we can import into the
    database as model info files.
    """
    pass


@model_info.command("crw")
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def crw_model_info(filename, output):
    """
    Parse the CRW metadata file and produce
    """
    r2dt.write_crw(filename, output)


@model_info.command("ribovision")
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def ribovision_model_info(filename, output):
    """
    Parse the metadata.tsv file from R2DT for Ribovision models to
    produce something we can put in our database.
    """
    r2dt.write_ribovision(filename, output)


@model_info.command("gtrnadb")
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def gtrnadb_model_info(filename, output):
    """
    Parse the metadata.tsv file from R2DT for gtrnadb models to
    produce something we can put in our database.
    """
    r2dt.write_gtrnadb(filename, output)


@model_info.command("rnase-p")
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def rnase_p_model_info(filename, output):
    """
    Parse the metadata.tsv file from R2DT for Ribovision models to
    produce something we can put in our database.
    """
    r2dt.write_rnase_p(filename, output)


@cli.command("create-attempted")
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def r2dt_create_attempted(filename, output):
    attempted.r2dt(filename, output)


@cli.command("publish")
@click.option("--suffix", default="")
@click.option("--allow-missing", is_flag=True, default=False)
@click.argument("model_info", type=click.File("r"))
@click.argument(
    "directory",
    type=click.Path(
        writable=False,
        dir_okay=True,
        file_okay=False,
    ),
)
@click.argument(
    "output",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
def r2dt_publish(model_info, directory, output, allow_missing, suffix=""):
    r2dt.publish(
        model_info, directory, output, allow_missing=allow_missing, suffix=suffix
    )


@cli.command("prepare-s3")
@click.option("--allow-missing", is_flag=True, default=False)
@click.argument("model_info", type=click.File("r"))
@click.argument(
    "directory",
    type=click.Path(
        writable=False,
        dir_okay=True,
        file_okay=False,
    ),
)
@click.argument(
    "output",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
@click.argument("file_list", type=click.Path())
def r2dt_prepare_s3(model_info, directory, output, file_list, allow_missing):
    file_list = Path(file_list)
    output = Path(output)
    r2dt.prepare_s3(
        model_info, directory, output, file_list, allow_missing=allow_missing
    )
