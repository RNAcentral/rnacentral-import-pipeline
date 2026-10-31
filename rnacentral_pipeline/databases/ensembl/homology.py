# -*- coding: utf-8 -*-

"""
Copyright [2009-2026] EMBL-European Bioinformatics Institute
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

import csv
import gzip
import hashlib
import typing as ty
from pathlib import Path

from rnacentral_pipeline.databases.ensembl import mrna


def canonical(gff_file: Path) -> ty.Dict[str, str]:
    """
    Map each protein coding gene to its Ensembl canonical transcript.
    """
    return {
        gene: transcript
        for transcript, (gene, tags) in mrna.protein_coding(gff_file).items()
        if "Ensembl_canonical" in tags
    }


def rows(
    homology_files: ty.Iterable[Path], gff_files: ty.Iterable[Path]
) -> ty.Iterable[ty.Tuple[str, str]]:
    """
    Turn Ensembl's gene-level homology pairs into (homology_group, transcript)
    rows for ensembl_compara, one group per pair of canonical mRNAs. A pair is
    only kept when both genes are in the given GFF3 files, and each species'
    file lists the same pair again from its own side.
    """
    transcripts = {}
    for gff_file in gff_files:
        transcripts.update(canonical(gff_file))

    seen = set()
    for path in homology_files:
        with gzip.open(path, "rt") as raw:
            for row in csv.DictReader(raw, delimiter="\t"):
                pair = tuple(
                    sorted((row["ref_gene_stable_id"], row["query_gene_stable_id"]))
                )
                if pair in seen or not all(gene in transcripts for gene in pair):
                    continue
                seen.add(pair)
                group = hashlib.sha256("|".join(pair).encode("utf-8")).hexdigest()
                for gene in pair:
                    yield (group, transcripts[gene])
