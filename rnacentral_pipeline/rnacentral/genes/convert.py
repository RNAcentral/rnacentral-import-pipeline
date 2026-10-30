# -*- coding: utf-8 -*-

"""
Copyright [2009-2025] EMBL-European Bioinformatics Institute
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

import re
import tempfile
from pathlib import Path

import gffutils
import polars as pl
from sqlitedict import SqliteDict

from rnacentral_pipeline.rnacentral.genes import extract as data

insdc_so_lookup = {
    "ncRNA": "SO:0000655",
    "lncRNA": "SO:0001877",
    "RNase_P_RNA": "SO:0000386",
    "SRP_RNA": "SO:0000590",
    "Y_RNA": "SO:0000405",
    "hammerhead_ribozyme": "SO:0000380",
    "miRNA": "SO:0000276",
    "misc_RNA": "SO:0000673",
    "pre_miRNA": "SO:0001244",
    "rRNA": "SO:0000252",
    "ribozyme": "SO:0000374",
    "sRNA": "SO:0000655",  ## this ine isn't an INSDC term, not sure what to do with it
    "scaRNA": "SO:0002095",
    "snoRNA": "SO:0000275",
    "tRNA": "SO:0000253",
    "telomerase_RNA": "SO:0000390",
    "vault_RNA": "SO:0000404",
    "snRNA": "SO:0000274",
    "piRNA": "SO:0001035",
    "precursor_RNA": "SO:0000655",
    "scRNA": "SO:0000013",
    "antisense_RNA": "SO:0000644",
    "other": "SO:0000655",
    "guide_RNA": "SO:0000602",
    "autocatalytically_spliced_intron": "SO:0000588",
    "RNase_MRP_RNA": "SO:0000385",
}


def get_assembly(path: Path) -> str:
    with path.open("r") as raw:
        for line in raw:
            if not line.startswith("#"):
                break
            if line.startswith("#!genome-version"):
                parts = line.split(" ", 1)
                return parts[1].strip()

    raise ValueError(f"Could not find assembly id in {path}")


## circRNAs and piRNAs are not clustered into genes; these are the `type`s the GFF export writes
EXCLUDED_RNA_TYPES = {"circRNA", "piRNA"}


def keep_valid_start(feature):
    """
    gffutils index transform: a GFF start below 1 is invalid and makes gffutils' bin calculation
    return a set, which SQLite cannot store (and gffutils then hides that error behind a
    misleading `str has no attribute decode`). Some exports contain such features (a start of 0 on the first
    bases of a sequence). Index it with start 1, and keep the value from the file in `raw_start` so the
    region name, which includes the coordinates, is still what the file says.
    """
    if feature.start < 1:
        feature.attributes["raw_start"] = [str(feature.start)]
        feature.start = 1
    return feature


def file_start(feature):
    """The start as written in the file."""
    return (
        int(feature.attributes["raw_start"][0])
        if "raw_start" in feature.attributes
        else feature.start
    )


def load_coordinates(
    path: Path,
    taxid: int,
    regions_data: Path | None = None,
    assembly_id: str | None = None,
    min_chromosome_transcripts: int = data.MIN_CHROMOSOME_TRANSCRIPTS,
) -> pl.DataFrame:
    """
    Load transcript coordinates from a GFF3 file and return a Polars DataFrame.
    Args:
        path (Path): Path to the GFF3 file.
        taxid (int): Taxonomy ID for the transcripts.
        regions_data (Path): regions of the taxid, used to find each region's assembly
        assembly_id (str): stamp every transcript with this assembly instead of
            joining to regions_data. Needed to replay old releases, whose regions
            are no longer in the database.
        min_chromosome_transcripts (int): drop chromosomes with fewer transcripts than this
    Returns:
        pl.DataFrame: DataFrame with transcript coordinates.
    """
    if assembly_id is None and regions_data is None:
        raise ValueError("either assembly_id or regions_data is required")
    transcript_data = []

    with tempfile.NamedTemporaryFile() as tmp:
        db = gffutils.create_db(str(path), tmp.name, transform=keep_valid_start)
        for ncrna_gene in db.features_of_type(("transcript")):
            region_start = file_start(ncrna_gene)
            region_stop = ncrna_gene.stop
            chromosome = ncrna_gene.chrom
            strand = ncrna_gene.strand
            urs = ncrna_gene.attributes["Name"][0]
            rna_insdc_type = ncrna_gene.attributes["type"][0]
            if rna_insdc_type in EXCLUDED_RNA_TYPES:
                continue
            rna_so_type = insdc_so_lookup.get(rna_insdc_type, "SO:0000655")
            if "_" in urs:
                urs_taxid = urs
            else:
                urs_taxid = f"{urs}_{taxid}"

            exons = []
            for exon in sorted(db.children(ncrna_gene.id), key=lambda f: f.start):
                exon_start = file_start(exon)
                exon_stop = exon.stop
                exon_dict = {
                    "chromosome": chromosome,
                    "urs_taxid": urs_taxid,
                    # "region_id": "" ## This has to be collected by querying later
                    "region_start": region_start,
                    "region_stop": region_stop,
                    "exon_start": exon_start,
                    "exon_stop": exon_stop,
                    "strand": 1 if strand == "+" else -1,
                    "so_type": rna_so_type,
                }
                exons.append(exon_dict)
            if not urs_taxid.startswith("URS"):
                raise ValueError(
                    f"{path}: transcript {ncrna_gene.id} is not a URS_taxid"
                )
            if len(exons) == 0:
                region_id = (
                    f"{urs_taxid}@{chromosome}/{region_start}-{region_stop}:{strand}"
                )
                transcript_dict = {
                    "chromosome": chromosome,
                    "region_name": region_id,
                    "urs_taxid": urs_taxid,
                    # "region_id": "" ## This has to be collected by querying later
                    "region_start": region_start,
                    "region_stop": region_stop,
                    "exon_start": region_start,
                    "exon_stop": region_stop,
                    "strand": 1 if strand == "+" else -1,
                    "so_type": rna_so_type,
                }
                transcript_data.append(transcript_dict)
            else:
                exon_regions = ",".join(
                    [f"{e['exon_start']}-{e['exon_stop']}" for e in exons]
                )
                # the GFF ID is <urs_taxid>.<n>, not the region name, so rebuild it as the database does
                region_id = f"{urs_taxid}@{chromosome}/{exon_regions}:{strand}"
                for e in exons:
                    e["region_name"] = region_id
                transcript_data.extend(exons)

    transcript_dataframe = pl.DataFrame(transcript_data)
    if assembly_id is not None:
        transcript_dataframe = transcript_dataframe.with_columns(
            pl.lit(assembly_id).alias("assembly_id")
        )
    else:
        ## A region_name can exist on several assemblies, so the regions must
        ## be for the assembly of this file or transcripts are duplicated
        transcript_dataframe = data.add_assembly_region_ids(
            transcript_dataframe, regions_data
        )

    transcript_dataframe = transcript_dataframe.group_by(
        ["assembly_id", "chromosome", "region_name"], maintain_order=True
    ).agg(
        pl.col("urs_taxid").first(),
        pl.col("region_start").first(),
        pl.col("region_stop").first(),
        pl.col("exon_start"),
        pl.col("exon_stop"),
        pl.col("strand").first(),
        pl.col("so_type").first(),
    )

    return data.drop_sparse_chromosomes(
        transcript_dataframe, min_chromosome_transcripts
    )


def gff_to_polars(
    path: Path,
    taxid: int,
    regions_data: Path | None = None,
    assembly_id: str | None = None,
    min_chromosome_transcripts: int = data.MIN_CHROMOSOME_TRANSCRIPTS,
) -> pl.DataFrame:
    """
    Convert a GFF3 file to a Polars DataFrame containing transcript coordinates.

    Args:
        path (Path): Path to the GFF3 file.
        taxid (int): Taxonomy ID for the transcripts.
        conn_str (str): Database connection string to fetch assembly IDs.

    Returns:
        pl.DataFrame: DataFrame with transcript coordinates.
    """
    return load_coordinates(
        path, taxid, regions_data, assembly_id, min_chromosome_transcripts
    )
