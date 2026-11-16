from itertools import product
from pathlib import Path

import click

from rnacentral_pipeline import schemas, writers
from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.gocams.search import find_gocam
from rnacentral_pipeline.output_format import format_option, is_parquet
from rnacentral_pipeline.writers import GoCamMembershipWriter

_PARQUET_SCHEMAS = {
    "go_cam_memberships": schemas.GO_CAM_ANNOTATION,
}


@click.group("gocams")
def cli():
    """
    Commands for GO-CAMS data.
    """


@cli.command("search")
@click.argument("filename", type=click.Path())
@click.argument(
    "output",
    default=".",
    type=click.Path(
        writable=True,
        dir_okay=True,
        file_okay=False,
    ),
)
@format_option
def search_gocams(filename, output):
    """
    Take a list of URS_taxids with known GO annotations and
    search them for Go-Cam membership. Write a parquet file containing flags and
    GO-CAM id
    """

    with open(filename, "r") as go_annotated_urs_taxids:
        urs_taxid_gen = (l.strip() for l in go_annotated_urs_taxids if l.strip())
        go_cams = (
            data.GoCamMembership(urs_taxid=u, go_cam_id=g)
            for u, g in find_gocam(
                urs_taxid_gen
            )  ## List of go_cam_ids matching this URS_taxid
        )

        if is_parquet():
            with writers.build_parquet(
                writers.GoCamMembershipWriter, Path(output), _PARQUET_SCHEMAS
            ) as writer:
                writer.write(go_cams)
            return
        with writers.build(writers.GoCamMembershipWriter, Path(output)) as writer:
            writer.write(go_cams)
