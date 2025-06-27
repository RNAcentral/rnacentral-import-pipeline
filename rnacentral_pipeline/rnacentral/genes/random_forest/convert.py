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

import tempfile
from pathlib import Path

import gffutils
import polars as pl
from sqlitedict import SqliteDict

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
                return parts[1]

    raise ValueError(f"Could not find assembly id in {path}")


def load_coordinates(path: Path, taxid: int) -> SqliteDict:
    # assembly_id = get_assembly(path)
    assembly_id = ".".join(path.stem.split(".")[1:])
    transcript_data = []

    with tempfile.NamedTemporaryFile() as tmp:
        db = gffutils.create_db(str(path), tmp.name)
        for ncrna_gene in db.features_of_type(("transcript")):
            region_start = ncrna_gene.start
            region_stop = ncrna_gene.stop
            exons = []
            for exon in db.children(ncrna_gene.id):
                urs = ncrna_gene.attributes["Name"][0]
                rna_insdc_type = ncrna_gene.attributes["type"][0]
                rna_so_type = insdc_so_lookup[rna_insdc_type]
                chromosome = ncrna_gene.chrom
                strand = ncrna_gene.strand

                exon_start = exon.start
                exon_stop = exon.stop

                if "_" in urs:
                    urs_taxid = urs
                else:
                    urs_taxid = f"{urs}_{taxid}"

                exon_dict = {
                    "assembly_id": assembly_id,
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
            if len(exons) == 0:
                region_id = (
                    f"{urs_taxid}@{chromosome}/{region_start}-{region_stop}:{strand}"
                )
                transcript_dict = {
                    "assembly_id": assembly_id,
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
                region_id = f"{urs_taxid}@{chromosome}/{exon_regions}:{strand}"
                for e in exons:
                    e["region_name"] = region_id
                transcript_data.extend(exons)

    transcript_dataframe = pl.DataFrame(transcript_data)
    print(transcript_dataframe.columns)
    transcript_dataframe = transcript_dataframe.group_by(
        ["assembly_id", "chromosome", "region_name"], maintain_order=True
    ).agg(
        pl.col("urs_taxid").first(),
        pl.col("region_start").first(),
        pl.col("region_stop").first(),
        # pl.col("region_id").first(),
        pl.col("exon_start"),
        pl.col("exon_stop"),
        pl.col("strand").first(),
        pl.col("so_type").first(),
    )
    return transcript_dataframe


def gff_to_polars(path: Path, taxid: int) -> pl.DataFrame:
    """
    Convert a GFF3 file to a Polars DataFrame containing transcript coordinates.

    Args:
        path (Path): Path to the GFF3 file.
        taxid (int): Taxonomy ID for the transcripts.

    Returns:
        pl.DataFrame: DataFrame with transcript coordinates.
    """
    return load_coordinates(path, taxid)
