# -*- coding: utf-8 -*-

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

import json
from pathlib import Path

import click
import polars as pl

from rnacentral_pipeline.rnacentral.genes import convert as gff_convert
from rnacentral_pipeline.rnacentral.genes import extract as data
from rnacentral_pipeline.rnacentral.genes import lineage, metadata, store


@click.group("genes")
def cli():
    """
    A group of commands dealing with building genes.
    """
    pass


@cli.group("utils")
def utils_cli():
    """
    Utility commands for gene processing.
    """
    pass


## Extract transcripts for the GNN gene model
@utils_cli.command("fetch")
@click.option(
    "--conn_str",
    envvar="PGDATABASE",
    required=True,
    help="Database connection string (can use PGDATABASE env var)",
)
@click.option("--taxid", type=int, required=True, help="Taxonomy ID to fetch data for")
@click.option(
    "--output_dir",
    required=True,
    type=click.Path(),
    help="Directory for transcripts__<assembly>.pq, one parquet per assembly",
)
@click.option(
    "--min_chromosome_transcripts",
    type=int,
    default=data.MIN_CHROMOSOME_TRANSCRIPTS,
    show_default=True,
    help="Drop chromosomes/scaffolds with fewer transcripts than this",
)
def fetch(conn_str, taxid, output_dir, min_chromosome_transcripts):
    """
    Fetch the transcripts of a taxon from the database, one parquet file per assembly.

    The GNN keys its nodes on region_name, which can be the same on two assemblies, so
    it must be run on one assembly at a time. Writes nothing for a taxon with no
    transcripts.
    """
    click.echo(f"Fetching data for taxid {taxid}...")
    transcripts = data.fetch_data(taxid, conn_str, min_chromosome_transcripts)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if transcripts.height == 0:
        click.echo("No transcripts")
        return
    for (assembly,), frame in transcripts.group_by("assembly_id"):
        safe = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in assembly)
        frame.write_parquet(out / f"transcripts__{safe}.pq")
        click.echo(f"Saved {frame.height} transcripts for {assembly}")


@utils_cli.command("convert")
@click.option("--gff_file", required=True, help="Input GFF file with genes")
@click.option("--taxid", type=int, help="Taxonomy ID for the GFF file")
@click.option(
    "--regions_data",
    default=None,
    help="Path to extracted region data for this taxid, ndjson. Not needed with --assembly_id",
)
@click.option(
    "--assembly_id",
    default=None,
    help="Assembly of this GFF. Skips the regions lookup, so old releases can be replayed",
)
@click.option(
    "--min_chromosome_transcripts",
    type=int,
    default=data.MIN_CHROMOSOME_TRANSCRIPTS,
    show_default=True,
    help="Drop chromosomes/scaffolds with fewer transcripts than this",
)
def convert(gff_file, taxid, regions_data, assembly_id, min_chromosome_transcripts):
    """
    Convert GFF file to parquet format.

    This command reads a GFF file containing gene annotations and converts it to a
    parquet format suitable for further processing.
    """
    if not Path(gff_file).exists():
        raise click.ClickException(f"GFF file not found: {gff_file}")
    if not regions_data and not assembly_id:
        raise click.ClickException("one of --regions_data or --assembly_id is required")
    output_name = Path(gff_file).stem + "_transcripts.parquet"
    output_path = Path(gff_file).parent / output_name
    transcripts_table = gff_convert.gff_to_polars(
        Path(gff_file),
        taxid=taxid,
        regions_data=Path(regions_data) if regions_data else None,
        assembly_id=assembly_id,
        min_chromosome_transcripts=min_chromosome_transcripts,
    )
    transcripts_table.write_parquet(output_path)
    click.echo(
        f"Converted {gff_file} to {output_path} with {transcripts_table.height} transcripts."
    )


@utils_cli.command("clusters")
@click.option("--transcripts", required=True, type=click.Path(exists=True))
@click.option("--assignments", required=True, type=click.Path(exists=True))
@click.option(
    "--assembly_id", default=None, help="Only use transcripts of this assembly"
)
@click.option("--output", required=True, type=click.Path())
def clusters(transcripts, assignments, assembly_id, output):
    """
    Turn the GNN's region_name -> gene_id assignments into a clusters parquet for `plan`.
    """
    result = lineage.clusters_from_assignments(
        pl.read_parquet(transcripts), pl.read_parquet(assignments), assembly_id
    )
    result.write_parquet(output)
    click.echo(f"Wrote {result.height} clustered regions to {output}")


@utils_cli.command("plan")
@click.option(
    "--state",
    default="-",
    help="State parquet from the previous release, or - for none",
)
@click.option(
    "--clusters",
    "cluster_files",
    multiple=True,
    type=click.Path(exists=True),
    help="Clusters of each assembly; none means the taxon has no genes (needs --allow_empty)",
)
@click.option("--release", type=int, required=True)
@click.option(
    "--prefix", default=None, help="Ensembl stable prefix; needed only to create genes"
)
@click.option(
    "--meta",
    default=None,
    type=click.Path(exists=True),
    help="meta.json from fetch-state; supplies the prefix",
)
@click.option("--out_dir", required=True, type=click.Path())
@click.option(
    "--allow_empty", is_flag=True, help="Allow no clusters at all (retires every gene)"
)
@click.option("--max_deactivate_fraction", type=float, default=0.5, show_default=True)
@click.option("--min_inherit", type=float, default=0.5, show_default=True)
@click.option("--min_relation", type=float, default=0.2, show_default=True)
@click.option(
    "--min_new_fraction",
    type=float,
    default=0.05,
    show_default=True,
    help="An old gene inherits only if it supplied this share of the cluster...",
)
@click.option(
    "--min_support",
    type=int,
    default=20,
    show_default=True,
    help="...or at least this many members",
)
def plan(
    state,
    cluster_files,
    release,
    prefix,
    meta,
    out_dir,
    allow_empty,
    max_deactivate_fraction,
    min_inherit,
    min_relation,
    min_new_fraction,
    min_support,
):
    """
    Work out how one release's genes follow on from the previous release's.

    Writes genes, state (feed it to the next release), events, member_deltas and
    stats.json to --out_dir. Nothing is written to the database.
    """
    previous = lineage.empty_state() if state == "-" else pl.read_parquet(state)
    if prefix is None and meta:
        prefix = json.loads(Path(meta).read_text())["prefix"]
    if cluster_files:
        found = pl.concat([pl.read_parquet(f) for f in cluster_files])
    else:
        found = pl.DataFrame(schema={c: pl.Utf8 for c in lineage.CLUSTER_COLUMNS})
    result = lineage.plan(
        previous,
        found,
        release,
        prefix,
        allow_empty=allow_empty,
        max_deactivate_fraction=max_deactivate_fraction,
        min_inherit=min_inherit,
        min_relation=min_relation,
        min_new_fraction=min_new_fraction,
        min_support=min_support,
    )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    result.genes.write_parquet(out / "genes.parquet")
    lineage.as_state(result.genes).write_parquet(out / "state.parquet")
    result.events.write_parquet(out / "events.parquet")
    result.member_deltas.write_parquet(out / "member_deltas.parquet")
    (out / "stats.json").write_text(json.dumps(result.stats))
    click.echo(f"Release {release}: {result.stats}")


@utils_cli.command("collect")
@click.option(
    "--plan_dir",
    "plan_dirs",
    multiple=True,
    required=True,
    type=click.Path(exists=True),
    help="Plan directories from `plan`, oldest release first",
)
@click.option("--output_dir", required=True, type=click.Path())
def collect(plan_dirs, output_dir):
    """
    Gather the plans of a chain of releases into one set of files for a bulk load.

    Writes events.parquet and member_deltas.parquet (all releases), and the genes and
    state of the last release. The result has the file names of a plan, so it can be
    loaded with `apply --include_inactive`.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name in ("events", "member_deltas"):
        pl.concat(
            [pl.read_parquet(Path(d) / f"{name}.parquet") for d in plan_dirs]
        ).write_parquet(out / f"{name}.parquet")
    last = Path(plan_dirs[-1])
    pl.read_parquet(last / "genes.parquet").write_parquet(out / "genes.parquet")
    pl.read_parquet(last / "state.parquet").write_parquet(out / "state.parquet")
    click.echo(f"Collected {len(plan_dirs)} releases into {out}")


@utils_cli.command("fetch-state")
@click.option("--taxid", type=int, required=True)
@click.option(
    "--output", required=True, type=click.Path(), help="State parquet for `plan`"
)
@click.option(
    "--meta", required=True, type=click.Path(), help="JSON with the Ensembl prefix"
)
@click.option("--db_str", envvar="PGDATABASE")
def fetch_state(taxid, output, meta, db_str):
    """
    Read the stored genes of a taxon, with their member history, as the state for `plan`.
    Read-only.
    """
    state, info = store.fetch_state(taxid, db_str)
    state.write_parquet(output)
    Path(meta).write_text(json.dumps(info))
    click.echo(f"Fetched {state.height} genes for taxid {taxid} ({info})")


@utils_cli.command("metadata")
@click.option(
    "--genes",
    required=True,
    type=click.Path(exists=True),
    help="genes.parquet from `plan`",
)
@click.option("--taxid", type=int, required=True)
@click.option("--output", required=True, type=click.Path())
@click.option("--db_str", envvar="PGDATABASE")
def compute_metadata(genes, taxid, output, db_str):
    """
    Pick a description and RNA type for every active gene. Reads evidence from the
    database in bulk and scores it here; writes nothing to the database.
    """
    active = pl.read_parquet(genes).filter(pl.col("is_active"))
    members = (
        active.select("internal_name", pl.col("members").alias("region_name"))
        .explode("region_name")
        .drop_nulls("region_name")
    )
    result, missing = metadata.compute_metadata(
        members, store.fetch_evidence(taxid, db_str)
    )
    result.write_parquet(output)
    if missing.height:
        click.echo(
            f"{missing.height} genes have no evidence, e.g. {missing['internal_name'][0]}"
        )
    click.echo(f"Wrote metadata for {result.height} genes to {output}")


@utils_cli.command("apply")
@click.option("--plan_dir", required=True, type=click.Path(exists=True))
@click.option("--metadata", "metadata_file", default=None, type=click.Path(exists=True))
@click.option("--taxid", type=int, required=True)
@click.option("--release", type=int, required=True)
@click.option(
    "--allow_unresolved",
    is_flag=True,
    help="Load even if some members are not in rnc_sequence_regions",
)
@click.option(
    "--include_inactive",
    is_flag=True,
    help="Also write genes that ended inactive; for loading a replayed history",
)
@click.option("--db_str", envvar="PGDATABASE")
def apply(
    plan_dir, metadata_file, taxid, release, allow_unresolved, include_inactive, db_str
):
    """
    Write a plan (genes, members, history, member history, metadata) in ONE transaction.
    """
    meta = pl.read_parquet(metadata_file) if metadata_file else None
    result = store.apply_plan(
        plan_dir, taxid, release, db_str, meta, allow_unresolved, include_inactive
    )
    click.echo(f"Applied plan for taxid {taxid}: {result}")
