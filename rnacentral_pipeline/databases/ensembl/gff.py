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

import tempfile
import typing as ty
from pathlib import Path

import gffutils
from sqlitedict import SqliteDict

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.ensembl.data import TranscriptInfo

SO_MAPPING = {
    "RNase_MRP_RNA": "SO:0000590",
    "RNase_P_RNA": "SO:0000386",
    "SRP_RNA": "SO:0000590",
    "Y_RNA": "SO:0000405",
    "guide_RNA": "SO:0000602",
    "lnc_RNA": "SO:0001877",
    "miRNA": "SO:0000276",
    "ncRNA": "SO:0000655",
    "ncRNA_gene": "SO:0001263",
    "pre_miRNA": "SO:0001244",
    "rRNA": "SO:0000252",
    "scRNA": "SO:0000013",
    "snRNA": "SO:0000274",
    "snoRNA": "SO:0000275",
    "tRNA": "SO:0000253",
    "telomerase_RNA": "SO:0000390",
    "tmRNA": "SO:0000584",
    "transcript": "SO:0000655",
    "processed_transcript": "SO:0001503",
}

IGNORED_TRANSCRIPTS = {
    "pseudogenic_transcript",
    "C_gene_segment",
    "D_gene_segment",
    "J_gene_segment",
    "V_gene_segment",
    "biological_region",
    "chromosome",
    "five_prime_UTR",
    "gene",
    "mRNA",
    "pseudogene",
    "pseudogenic_transcript",
    "scaffold",
    "three_prime_UTR",
    "unconfirmed_transcript",
}


def get_assembly(path: Path) -> str:
    path = Path(path)
    with path.open("r") as raw:
        for line in raw:
            if not line.startswith("#"):
                break
            if line.startswith("#!genome-version"):
                parts = line.split(" ", 1)
                return parts[1].strip()

    raise ValueError(f"Could not find assembly id in {path}")


def get_assembly_accession(path: Path) -> str:
    with path.open("r") as raw:
        for line in raw:
            if not line.startswith("#"):
                break
            if line.startswith("#!genome-build-accession"):
                return line.split(" ", 1)[1].strip()

    raise ValueError(f"Could not find assembly accession in {path}")


def load_coordinates(path: Path, assembly_id: ty.Optional[str] = None) -> SqliteDict:
    assembly_id = assembly_id or get_assembly(path)
    mapping = SqliteDict()
    # dir="." because singularity --contain gives us a tiny in-memory /tmp
    with tempfile.NamedTemporaryFile(dir=".") as tmp:
        db = gffutils.create_db(str(path), tmp.name)
        for ncrna_gene in db.features_of_type("ncRNA_gene"):
            transcripts = [
                t
                for t in db.children(ncrna_gene)
                if t.featuretype != "exon" and t.featuretype not in IGNORED_TRANSCRIPTS
            ]
            # Some genesets list an ncRNA as a gene alone, with no transcript or exon.
            if not transcripts:
                pid = ncrna_gene["ID"][0].split(":", 1)[1]
                biotype = ncrna_gene.attributes.get("biotype", ["ncRNA"])[0]
                mapping[pid] = TranscriptInfo(
                    so_rna_type=SO_MAPPING.get(biotype, SO_MAPPING["ncRNA"]),
                    regions=[
                        data.SequenceRegion(
                            assembly_id=assembly_id,
                            chromosome=ncrna_gene.chrom,
                            strand=ncrna_gene.strand,
                            exons=[
                                data.Exon(start=ncrna_gene.start, stop=ncrna_gene.stop)
                            ],
                            coordinate_system=data.CoordinateSystem.one_based(),
                        )
                    ],
                    from_gencode=False,
                )
            for transcript in transcripts:
                exons: ty.List[data.Exon] = []
                for exon in db.children(
                    transcript, featuretype="exon", order_by="start"
                ):
                    exons.append(data.Exon(start=exon.start, stop=exon.stop))

                gencode = "havana" in transcript.source
                pid = transcript["ID"][0].split(":", 1)[1]
                if pid in mapping:
                    raise ValueError(f"Duplicate ensembl ID seen {pid}")

                mapping[pid] = TranscriptInfo(
                    so_rna_type=SO_MAPPING[transcript.featuretype],
                    regions=[
                        data.SequenceRegion(
                            assembly_id=assembly_id,
                            chromosome=ncrna_gene.chrom,
                            strand=transcript.strand,
                            exons=exons,
                            coordinate_system=data.CoordinateSystem.one_based(),
                        )
                    ],
                    from_gencode=gencode,
                )
    mapping.commit()
    return mapping
