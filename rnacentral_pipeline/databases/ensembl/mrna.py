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

import re
import typing as ty
from pathlib import Path

import attr
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.ensembl import gff
from rnacentral_pipeline.databases.ensembl import helpers as common
from rnacentral_pipeline.databases.ensembl.parser import URL
from rnacentral_pipeline.databases.ensembl.vertebrates import helpers
from rnacentral_pipeline.databases.helpers import embl

MRNA = "SO:0000234"
FIVE_PRIME_UTR = "SO:0000204"
THREE_PRIME_UTR = "SO:0000205"

TAGS = {"Ensembl_canonical", "MANE_Select", "MANE_Plus_Clinical"}

XREFS = {
    "CCDS",
    "RefSeq_mRNA",
    "RefSeq_mRNA_predicted",
    "RefSeq_peptide",
    "RefSeq_peptide_predicted",
    "Uniprot/SPTREMBL",
    "Uniprot/SWISSPROT",
    "Uniprot_isoform",
}

UTRS = {
    FIVE_PRIME_UTR: ("five_prime_UTR", "five_prime_utr", "5' UTR"),
    THREE_PRIME_UTR: ("three_prime_UTR", "three_prime_utr", "3' UTR"),
}


def protein_coding(gff_file: Path) -> ty.Dict[str, ty.Tuple[str, ty.List[str]]]:
    """
    Map each protein_coding transcript to its gene and kept tags. The EMBL mRNA
    features carry no biotype, and also cover NMD and non-stop decay
    transcripts, so these come from the GFF3.
    """
    transcripts = {}
    with gff_file.open("r") as raw:
        for line in raw:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue
            attrs = dict(a.split("=", 1) for a in parts[8].split(";") if "=" in a)
            if attrs.get("biotype") == "protein_coding":
                tags = [t for t in attrs.get("tag", "").split(",") if t in TAGS]
                gene = attrs["Parent"].split(":", 1)[1]
                transcripts[attrs["ID"].split(":", 1)[1]] = (gene, sorted(tags))
    return transcripts


def xref_data(gene, cds) -> ty.Dict[str, ty.List[str]]:
    xrefs = {k: v for k, v in embl.xref_data(cds).items() if k in XREFS}
    # HGNC is only in the gene note, e.g. [Source:HGNC Symbol;Acc:HGNC:29173]
    hgnc = re.search(r"Acc:(HGNC:\d+)", " ".join(helpers.notes(gene)))
    if hgnc:
        xrefs["HGNC"] = [hgnc.group(1)]
    return xrefs


def clip(parts, start: int, end: int):
    kept = [
        FeatureLocation(max(int(p.start), start), min(int(p.end), end), p.strand)
        for p in parts
        if int(p.start) < end and int(p.end) > start
    ]
    if not kept:
        return None
    if len(kept) == 1:
        return kept[0]
    return CompoundLocation(kept)


def utr_locations(record, mrna, cds):
    """
    The UTRs are the parts of the mRNA exons outside the CDS span. Parts stay
    in transcript order, so extracting them gives the spliced sequence.
    """
    parts = mrna.location.parts
    cds_start, cds_end = int(cds.location.start), int(cds.location.end)
    upstream = clip(parts, 0, cds_start)
    downstream = clip(parts, cds_end, len(record))
    if parts[0].strand == -1:
        return {FIVE_PRIME_UTR: downstream, THREE_PRIME_UTR: upstream}
    return {FIVE_PRIME_UTR: upstream, THREE_PRIME_UTR: downstream}


def as_entry(record, gene, feature, accession, rna_type, label, gca) -> data.Entry:
    species, common_name = helpers.organism_naming(record)
    name = embl.locus_tag(gene) or embl.gene(gene).split(".")[0]
    organism = f"{species} ({common_name})" if common_name else species
    transcript = helpers.accession(feature)
    return data.Entry(
        primary_id=helpers.primary_id(feature),
        accession=accession,
        ncbi_tax_id=embl.taxid(record),
        database="ENSEMBL_MRNA",
        sequence=embl.sequence(record, feature),
        regions=common.regions(record, feature),
        rna_type=rna_type,
        url=URL.format(accession=gca, transcript=helpers.primary_id(feature)),
        seq_version=transcript.split(".", 1)[1] if "." in transcript else "1",
        lineage=embl.lineage(record),
        chromosome=helpers.chromosome(record),
        parent_accession=record.id,
        common_name=common_name,
        species=species,
        gene=embl.gene(gene),
        locus_tag=embl.locus_tag(gene),
        optional_id=embl.gene(gene),
        references=helpers.references(),
        mol_type="genomic DNA",
        description=f"{organism} {name} {label}",
    )


def transcript_entries(record, gene, mrna, cds, tags, gca) -> ty.Iterable[data.Entry]:
    accession = helpers.accession(mrna)
    entry = as_entry(record, gene, mrna, accession, MRNA, "mRNA", gca)
    # 1-based inclusive genomic span, so the export can draw the CDS exactly
    # even where a UTR was too short to import
    note = {"cds": [int(cds.location.start) + 1, int(cds.location.end)]}
    codon_start = int(cds.qualifiers.get("codon_start", ["1"])[0])
    if codon_start != 1:
        note["cds_phase"] = codon_start - 1
    if tags:
        note["tags"] = tags
    entry = attr.evolve(entry, xref_data=xref_data(gene, cds), note_data=note)
    if not entry.is_valid():
        return

    utrs = []
    for rna_type, location in utr_locations(record, mrna, cds).items():
        if location is None:
            continue
        suffix, _, label = UTRS[rna_type]
        feature = SeqFeature(location, type=suffix, qualifiers=mrna.qualifiers)
        utr = as_entry(
            record, gene, feature, f"{accession}:{suffix}", rna_type, label, gca
        )
        if utr.is_valid():
            utrs.append(utr)

    related = [
        data.RelatedSequence(sequence_id=u.accession, relationship=UTRS[u.rna_type][1])
        for u in utrs
    ]
    yield attr.evolve(entry, related_sequences=related)
    for utr in utrs:
        back = data.RelatedSequence(sequence_id=accession, relationship="mrna")
        yield attr.evolve(utr, related_sequences=[back])


def parse(raw: ty.IO, gff_file: Path) -> ty.Iterable[data.Entry]:
    gca = gff.get_assembly_accession(gff_file)
    wanted = protein_coding(gff_file)
    for record in SeqIO.parse(raw, "embl"):
        gene = None
        mrnas = {}
        for feature in record.features:
            if embl.is_gene(feature):
                gene = feature
            elif feature.type == "mRNA":
                mrnas[helpers.primary_id(feature)] = (gene, feature)
            elif feature.type == "CDS":
                pid = helpers.primary_id(feature)
                if pid in wanted and pid in mrnas:
                    mrna_gene, mrna = mrnas.pop(pid)
                    tags = wanted[pid][1]
                    yield from transcript_entries(
                        record, mrna_gene, mrna, feature, tags, gca
                    )
