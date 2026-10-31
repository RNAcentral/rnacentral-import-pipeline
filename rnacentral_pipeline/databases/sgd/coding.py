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

import collections
import typing as ty
from urllib.parse import unquote

import attr
from Bio import SeqIO
from Bio.Seq import Seq

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.ensembl.mrna import (
    FIVE_PRIME_UTR,
    MRNA,
    THREE_PRIME_UTR,
    UTRS,
)
from rnacentral_pipeline.databases.helpers import publications as pubs

TAXID = 559292
# SGD's genome has not changed since R64-1-1, the assembly RNAcentral uses
ASSEMBLY = "R64-1-1"
URL = "https://www.yeastgenome.org/locus/{sgdid}"
INTRONS = {"intron", "five_prime_UTR_intron"}
# The SGD database paper, and Pelechano 2013 for the transcript ends
REFERENCES = [22110037, 23615609]


def attributes(column: str) -> ty.Dict[str, str]:
    pairs = (a.split("=", 1) for a in column.split(";") if "=" in a)
    return {key: unquote(value) for key, value in pairs}


def chromosome(name: str) -> str:
    name = name[3:] if name.startswith("chr") else name
    return "Mito" if name == "mt" else name


def exons(start: int, stop: int, introns) -> ty.List[ty.Tuple[int, int]]:
    result = []
    for intron_start, intron_stop in sorted(introns):
        result.append((start, intron_start - 1))
        start = intron_stop + 1
    result.append((start, stop))
    # Without a measured transcript end, the mRNA can stop at a 5' UTR intron
    return [(s, e) for s, e in result if s <= e]


def clip(parts, start: int, stop: int) -> ty.List[ty.Tuple[int, int]]:
    return [(max(s, start), min(e, stop)) for s, e in parts if s <= stop and e >= start]


def as_entry(genome, mrna, gene, parts, accession, rna_type, label) -> data.Entry:
    chrom, strand = mrna["chromosome"], mrna["strand"]
    sequence = "".join(str(genome[chrom][s - 1 : e]) for s, e in parts).upper()
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    name = gene.get("gene") or gene["ID"]
    return data.Entry(
        primary_id=accession,
        accession=f"SGD_CODING:{accession}",
        ncbi_tax_id=TAXID,
        database="SGD_CODING",
        sequence=sequence,
        regions=[
            data.SequenceRegion(
                assembly_id=ASSEMBLY,
                chromosome=chromosome(chrom),
                strand=strand,
                exons=[data.Exon(start=s, stop=e) for s, e in parts],
                coordinate_system=data.CoordinateSystem.one_based(),
            )
        ],
        rna_type=rna_type,
        url=URL.format(sgdid=gene["dbxref"].split(":", 1)[1]),
        seq_version="1",
        chromosome=chromosome(chrom),
        gene=gene.get("gene"),
        locus_tag=gene["ID"],
        optional_id=gene.get("gene"),
        references=[pubs.reference(pmid) for pmid in REFERENCES],
        mol_type="genomic DNA",
        description=f"Saccharomyces cerevisiae S288C {name} {label}",
    )


def transcript_entries(genome, mrna, gene, children) -> ty.Iterable[data.Entry]:
    cds = [(s, e) for kind, s, e in children if kind == "CDS"]
    if not cds:
        return
    cds_start, cds_stop = min(s for s, _ in cds), max(e for _, e in cds)
    if cds_start < mrna["start"] or cds_stop > mrna["stop"]:
        return
    introns = [(s, e) for kind, s, e in children if kind in INTRONS]
    parts = exons(mrna["start"], mrna["stop"], introns)
    accession = mrna["ID"]
    entry = as_entry(genome, mrna, gene, parts, accession, MRNA, "mRNA")
    if not entry.is_valid():
        return

    upstream = clip(parts, mrna["start"], cds_start - 1)
    downstream = clip(parts, cds_stop + 1, mrna["stop"])
    if mrna["strand"] == "-":
        upstream, downstream = downstream, upstream
    utrs = []
    for rna_type, utr_parts in (
        (FIVE_PRIME_UTR, upstream),
        (THREE_PRIME_UTR, downstream),
    ):
        if not utr_parts:
            continue
        suffix, _, label = UTRS[rna_type]
        utr = as_entry(
            genome, mrna, gene, utr_parts, f"{accession}:{suffix}", rna_type, label
        )
        if utr.is_valid():
            utrs.append(utr)

    related = [
        data.RelatedSequence(sequence_id=u.accession, relationship=UTRS[u.rna_type][1])
        for u in utrs
    ]
    yield attr.evolve(entry, related_sequences=related)
    for utr in utrs:
        back = data.RelatedSequence(sequence_id=entry.accession, relationship="mrna")
        yield attr.evolve(utr, related_sequences=[back])


def parse(handle: ty.IO) -> ty.Iterable[data.Entry]:
    genes = {}
    mrnas = {}
    children = collections.defaultdict(list)
    for line in handle:
        if line.startswith("##FASTA"):
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 9:
            continue
        attrs = attributes(parts[8])
        start, stop = int(parts[3]), int(parts[4])
        if parts[2] == "gene":
            genes[attrs["ID"]] = attrs
        elif parts[2] == "mRNA":
            mrnas[attrs["ID"]] = dict(
                attrs, chromosome=parts[0], start=start, stop=stop, strand=parts[6]
            )
        for parent in attrs.get("Parent", "").split(","):
            if parent:
                children[parent].append((parts[2], start, stop))
    genome = {record.id: record.seq for record in SeqIO.parse(handle, "fasta")}

    for mrna in mrnas.values():
        gene = genes[mrna["Parent"]]
        # SGD judges these unlikely to encode a protein
        if gene.get("orf_classification") == "Dubious":
            continue
        yield from transcript_entries(genome, mrna, gene, children[mrna["ID"]])
