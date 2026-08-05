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
from rnacentral_pipeline.rnacentral.r2dt import s3 as r2dt_s3


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
@click.argument("info", type=click.File("r"))
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def crw_model_info(info, filename, output):
    """
    Parse the CRW metadata file and produce
    """
    r2dt.write_crw(filename, output, extra=info)


@model_info.command("ribovision")
@click.argument("info", type=click.File("r"))
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def ribovision_model_info(info, filename, output):
    """
    Parse the metadata.tsv file from R2DT for Ribovision models to
    produce something we can put in our database.
    """
    r2dt.write_ribovision(filename, output, extra=info)


@model_info.command("gtrnadb")
@click.argument("filename", type=click.File())
@click.argument("output", default="-", type=click.File("w"))
def gtrnadb_model_info(filename, output):
    """
    Parse the metadata.tsv file from R2DT for gtrnadb models to
    produce something we can put in our database.
    """
    r2dt.write_gtrnadb(filename, output)


@model_info.command("rnase-p")
@click.argument("info", type=click.File("r"))
@click.argument("filename", type=click.File("r"))
@click.argument("output", default="-", type=click.File("w"))
def rnase_p_model_info(info, filename, output):
    """
    Parse the metadata.tsv file from R2DT for RNase P models to
    produce something we can put in our database.
    """
    r2dt.write_rnase_p(filename, output, extra=info)


@model_info.command("rfam")
@click.argument("filename", type=click.File("r"))
@click.argument("db_url", type=str)
@click.argument("output", default="-", type=click.File("w"))
def rnase_p_model_info(filename, db_url, output):
    """
    Parse the metadata.tsv file from R2DT for Ribovision models to
    produce something we can put in our database.
    """
    r2dt.write_rfam(filename, db_url, output)


@cli.command("create-attempted")
@click.argument("filename", type=click.File("r"))
@click.argument("version", type=click.File("r"))
@click.argument("output", type=click.Path())
def r2dt_create_attempted(filename, version, output):
    version_string = version.read().strip()
    attempted.r2dt(filename, version_string, output)


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


@cli.command("upload-s3")
@click.option("--env", default="prod", help="top-level prefix / environment")
@click.option("--workers", default=32, help="parallel uploads")
@click.option("--endpoint", default=r2dt_s3.ENDPOINT)
@click.option("--bucket", default=r2dt_s3.BUCKET)
@click.option(
    "--allow-failures",
    default=0,
    help="tolerate up to this many failed objects instead of failing the batch",
)
@click.option(
    "--failure-list",
    default=None,
    type=click.Path(dir_okay=False),
    help="write the URS of each failed upload here, for later reconciliation",
)
@click.argument("file_list", type=click.Path(exists=True, dir_okay=False))
def r2dt_upload_s3(
    file_list, env, workers, endpoint, bucket, allow_failures, failure_list
):
    """
    Upload the gzipped SVGs listed in FILE_LIST (produced by prepare-s3) to S3.

    Replaces update-svg.sh: uploads via boto3 with server-side Content-MD5
    verification and fails loudly if any object does not upload, rather than
    silently swallowing HTTP errors.

    --allow-failures lets a stray unwritable object through rather than losing
    the batch; anything beyond the cap still fails. Use --failure-list to record
    what was skipped.
    """
    r2dt_s3.upload(
        file_list,
        env,
        endpoint=endpoint,
        bucket=bucket,
        workers=workers,
        allow_failures=allow_failures,
        failure_list=failure_list,
    )


@cli.command("drop-failed-uploads")
@click.option(
    "--max-failures",
    default=None,
    type=int,
    help="fail if FAILURE_LIST holds more URS than this (a run-wide cap)",
)
@click.argument("failure_list", type=click.Path(dir_okay=False))
@click.argument("data_files", nargs=-1, type=click.Path(exists=True, dir_okay=False))
def r2dt_drop_failed_uploads(failure_list, data_files, max_failures):
    """
    Remove rows for URS listed in FAILURE_LIST from DATA_FILES, in place.

    Keeps the database honest when upload-s3 tolerated a failure: the URS would
    otherwise get a structure row pointing at an object that is not in S3. Pass
    the hits and the attempted files both, so the sequence is left undone rather
    than recorded as attempted and never retried.

    A missing or empty FAILURE_LIST is a no-op, so callers can pass it
    unconditionally. --max-failures caps the whole run, since upload-s3's own
    tolerance is per batch.
    """
    r2dt_s3.drop_failed(failure_list, data_files, max_failures=max_failures)


@cli.command("verify-s3")
@click.option("--env", default="prod")
@click.option("--workers", default=32)
@click.option("--endpoint", default=r2dt_s3.ENDPOINT)
@click.option("--bucket", default=r2dt_s3.BUCKET)
@click.argument("file_list", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", default="-", type=click.File("w"))
def r2dt_verify_s3(file_list, output, env, workers, endpoint, bucket):
    """
    Checksum the files in FILE_LIST against what is actually in the bucket.

    Writes one line per problem (MISSING / SIZE / CHECKSUM / ERROR) to OUTPUT and
    exits non-zero if anything is missing or does not match.
    """
    problems = r2dt_s3.verify(
        file_list, env, endpoint=endpoint, bucket=bucket, workers=workers, out=output
    )
    if problems:
        raise click.ClickException(
            f"{problems} files missing or mismatched in {bucket}"
        )


@cli.command("list-s3")
@click.option("--env", default="prod")
@click.option("--workers", default=32)
@click.option("--depth", default=2, help="delimiter levels to shard on")
@click.option("--endpoint", default=r2dt_s3.ENDPOINT)
@click.option("--bucket", default=r2dt_s3.BUCKET)
@click.argument("output", default="-", type=click.File("w"))
def r2dt_list_s3(output, env, workers, depth, endpoint, bucket):
    """
    List the bare URS of every SVG in the bucket (one per line) to OUTPUT.

    Used to reconcile what is in S3 against the DB should_show set.
    """
    r2dt_s3.list_svgs(
        env, output, endpoint=endpoint, bucket=bucket, workers=workers, depth=depth
    )


@cli.command("sync-s3")
@click.option("--env", default="prod")
@click.option("--db-url", envvar="PGDATABASE")
@click.option("--workers", default=32)
@click.option("--depth", default=2, help="delimiter levels to shard on")
@click.option("--endpoint", default=r2dt_s3.ENDPOINT)
@click.option("--bucket", default=r2dt_s3.BUCKET)
@click.option("--missing-list", default="missing-svgs.txt")
@click.option("--orphan-list", default="orphan-svgs.txt")
@click.option(
    "--delete",
    is_flag=True,
    default=False,
    help="actually remove the orphans; without it this is a dry run",
)
def r2dt_sync_s3(
    env, db_url, workers, depth, endpoint, bucket, missing_list, orphan_list, delete
):
    """
    Reconcile the bucket against the database should-show set.

    Writes MISSING-LIST (should show, but no SVG in S3 -- feed this back through
    r2dt) and ORPHAN-LIST (in S3, but no longer should show). Nothing is removed
    unless --delete is given: check the orphan list first, since a mis-set --env
    would prune the wrong environment.
    """
    counts = r2dt_s3.sync(
        env,
        db_url,
        missing_list,
        orphan_list,
        delete=delete,
        endpoint=endpoint,
        bucket=bucket,
        workers=workers,
        depth=depth,
    )
    click.echo(
        f"in_s3={counts['in_s3']} missing={counts['missing']} "
        f"orphan={counts['orphan']}"
    )


@cli.command("download-s3")
@click.option("--env", default="prod")
@click.option("--workers", default=32)
@click.option("--depth", default=2, help="delimiter levels to shard on")
@click.option("--endpoint", default=r2dt_s3.ENDPOINT)
@click.option("--bucket", default=r2dt_s3.BUCKET)
@click.option(
    "--urs-list",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="URS to fetch, one per line; defaults to everything in the bucket",
)
@click.argument("directory", type=click.Path(file_okay=False))
def r2dt_download_s3(directory, env, workers, depth, endpoint, bucket, urs_list):
    """
    Download every SVG into DIRECTORY, ready to tar + gzip for distribution.

    Files are written decompressed as URS/xx/xx/xx/xx/<urs>.svg, mirroring the S3
    layout, with a manifest.tsv of urs/path/size/md5 alongside. Re-running skips
    files already on disk, so an interrupted run can just be repeated; the
    manifest is appended to, so it may then hold duplicate rows.
    """
    r2dt_s3.download(
        directory,
        env,
        urs_list=urs_list,
        endpoint=endpoint,
        bucket=bucket,
        workers=workers,
        depth=depth,
    )


@cli.command("prepare-sequences")
@click.argument("xref_urs", type=click.File("r"))
@click.argument("tracked_urs", type=click.File("r"))
@click.argument("urs_to_fetch", type=click.File("w"))
@click.option("--max_sequences", default=-1)
def r2dt_prepare_sequences(xref_urs, tracked_urs, urs_to_fetch, max_sequences):
    """
    Prepare a list of URS identifiers to fetch sequences for.

    This takes a file of all URS identifiers from cross-references and a file
    of already tracked URS identifiers. It produces a file of URS identifiers
    that are in the xref file but not in the tracked file. This can be limited
    to a maximum number of sequences.
    """
    r2dt.prepare_sequences(xref_urs, tracked_urs, urs_to_fetch, max_sequences)
