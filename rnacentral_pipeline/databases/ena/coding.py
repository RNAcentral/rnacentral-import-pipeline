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

import typing as ty

import attr
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.ensembl.mrna import (
    FIVE_PRIME_UTR,
    MRNA,
    THREE_PRIME_UTR,
    UTRS,
)
from rnacentral_pipeline.databases.helpers import embl

URL = "https://www.ebi.ac.uk/ena/browser/view/{accession}"

Parts = ty.List[ty.Tuple[int, int]]


def chromosome(record) -> str:
    source = record.features[0].qualifiers
    if "chromosome" in source:
        return source["chromosome"][0]
    # Ensembl's names for a bacterial chromosome and its plasmids
    if "plasmid" in source:
        return f"{source['plasmid'][0]}_plasmid"
    return "Chromosome"


def parts(loc) -> Parts:
    return sorted((int(p.start) + 1, int(p.end)) for p in loc.parts)


def clip(spans: Parts, start: int, stop: int) -> Parts:
    return [(max(s, start), min(e, stop)) for s, e in spans if s <= stop and e >= start]


def location(spans: Parts, strand: int):
    locs = [FeatureLocation(s - 1, e, strand) for s, e in spans]
    if strand == -1:
        locs.reverse()
    return locs[0] if len(locs) == 1 else CompoundLocation(locs)


def first(feature, name: str) -> ty.Optional[str]:
    return feature.qualifiers.get(name, [None])[0]


@attr.s(frozen=True)
class Transcript:
    accession: str = attr.ib()
    spans: Parts = attr.ib()
    strand: int = attr.ib()
    cds: ty.Tuple[int, int] = attr.ib()
    name: str = attr.ib()
    gene: ty.Optional[str] = attr.ib()
    locus_tag: ty.Optional[str] = attr.ib()
    url: str = attr.ib()


def whole_record(record) -> ty.Iterable[Transcript]:
    # Cut around the sense ORFs, as the probing references cut their UTRs
    cds = [f for f in record.features if f.type == "CDS" and f.location.strand == 1]
    if not cds:
        return
    starts = [int(f.location.start) + 1 for f in cds]
    stops = [int(f.location.end) for f in cds]
    name = record.description.split(",")[0].rstrip(".")
    yield Transcript(
        accession=record.id,
        spans=[(1, len(record))],
        strand=1,
        cds=(min(starts), max(stops)),
        name=name,
        gene=None,
        locus_tag=None,
        url=URL.format(accession=record.id),
    )


def annotated(record) -> ty.Iterable[Transcript]:
    cds = {}
    for feature in record.features:
        if feature.type == "CDS" and "pseudo" not in feature.qualifiers:
            cds[first(feature, "locus_tag")] = feature
    organism = record.annotations["organism"]

    mrnas = [f for f in record.features if f.type == "mRNA"]
    # Bacterial annotation has no mRNA features, so each CDS stands in for one
    if mrnas:
        pairs = [(m, cds.get(first(m, "locus_tag"))) for m in mrnas]
    else:
        pairs = [(c, c) for c in cds.values()]

    for feature, coding in pairs:
        if coding is None or "pseudo" in feature.qualifiers:
            continue
        coding_spans = parts(coding.location)
        locus_tag = first(feature, "locus_tag")
        gene = first(feature, "gene")
        transcript_id = first(feature, "transcript_id")
        accession = transcript_id or locus_tag
        url = URL.format(accession=first(coding, "protein_id") or record.id)
        yield Transcript(
            accession=accession,
            spans=parts(feature.location),
            strand=feature.location.strand,
            cds=(coding_spans[0][0], coding_spans[-1][1]),
            name=f"{organism} {gene or locus_tag}",
            gene=gene,
            locus_tag=locus_tag,
            url=url,
        )


def as_entry(
    record,
    transcript,
    spans,
    accession,
    rna_type,
    label,
    taxid,
    assembly,
    database,
    references,
) -> data.Entry:
    sequence = location(spans, transcript.strand).extract(record.seq)
    regions = []
    if assembly:
        regions = [
            data.SequenceRegion(
                assembly_id=assembly,
                chromosome=chromosome(record),
                strand="+" if transcript.strand == 1 else "-",
                exons=[data.Exon(start=s, stop=e) for s, e in spans],
                coordinate_system=data.CoordinateSystem.one_based(),
            )
        ]
    version = transcript.accession.rsplit(".", 1)
    return data.Entry(
        primary_id=accession,
        accession=f"{database}:{accession}",
        ncbi_tax_id=taxid,
        database=database,
        sequence=str(sequence).upper(),
        regions=regions,
        rna_type=rna_type,
        url=transcript.url,
        seq_version=version[1] if len(version) == 2 else "1",
        parent_accession=record.id,
        chromosome=chromosome(record) if assembly else None,
        gene=transcript.gene,
        locus_tag=transcript.locus_tag,
        optional_id=transcript.locus_tag,
        references=references,
        mol_type=first(record.features[0], "mol_type"),
        description=f"{transcript.name} {label}",
    )


def entries(
    record, transcript, taxid, assembly, database, references
) -> ty.Iterable[data.Entry]:
    start, stop = transcript.spans[0][0], transcript.spans[-1][1]
    cds_start, cds_stop = transcript.cds
    if cds_start < start or cds_stop > stop:
        return
    accession = transcript.accession
    entry = as_entry(
        record,
        transcript,
        transcript.spans,
        accession,
        MRNA,
        "mRNA",
        taxid,
        assembly,
        database,
        references,
    )
    # The IGV export draws the CDS from this 1-based genomic span, as for Ensembl
    entry = attr.evolve(entry, note_data={"cds": [cds_start, cds_stop]})
    if not entry.is_valid():
        return

    upstream = clip(transcript.spans, start, cds_start - 1)
    downstream = clip(transcript.spans, cds_stop + 1, stop)
    if transcript.strand == -1:
        upstream, downstream = downstream, upstream
    utrs = []
    for rna_type, utr_spans in (
        (FIVE_PRIME_UTR, upstream),
        (THREE_PRIME_UTR, downstream),
    ):
        if not utr_spans:
            continue
        suffix, _, label = UTRS[rna_type]
        utr = as_entry(
            record,
            transcript,
            utr_spans,
            f"{accession}:{suffix}",
            rna_type,
            label,
            taxid,
            assembly,
            database,
            references,
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


def cut(record, end: int):
    # Slicing drops the record's annotations and its whole-record source feature
    trimmed = record[:end]
    trimmed.annotations = record.annotations
    source = record.features[0]
    trimmed.features.insert(
        0,
        SeqFeature(
            FeatureLocation(0, end), type=source.type, qualifiers=source.qualifiers
        ),
    )
    return trimmed


def parse(
    handle: ty.IO,
    taxid: int,
    assembly: ty.Optional[str] = None,
    end: ty.Optional[int] = None,
) -> ty.Iterable[data.Entry]:
    # With an assembly the records are a genome; without one, virus genomes or segments
    for record in SeqIO.parse(handle, "embl"):
        if end:
            record = cut(record, end)
        # Ensembl has no organelle sequences for these genomes, so there is nowhere to place them
        if "organelle" in record.features[0].qualifiers:
            continue
        transcripts = annotated(record) if assembly else whole_record(record)
        references = embl.references(record)
        for transcript in transcripts:
            yield from entries(
                record, transcript, taxid, assembly, "ENA_CODING", references
            )
