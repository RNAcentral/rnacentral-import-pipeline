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
import json
import typing as ty

import attr
from Bio import SeqIO

from rnacentral_pipeline.databases import data
from rnacentral_pipeline.databases.data.utils import SO_INSDC_MAPPING, as_so_term
from rnacentral_pipeline.databases.ena import coding
from rnacentral_pipeline.databases.helpers import embl
from rnacentral_pipeline.databases.helpers import publications as pubs

URL = "https://regulondb.ccg.unam.mx/operon/{operon}"

# Ensembl's name for the U00096.3 genome that RegulonDB annotates
ASSEMBLY = "ASM584v2"

# Strong and Confirmed only: Weak and unrated TUs rest on indirect evidence
CONFIDENCE = {"S", "C"}

PRECURSOR_RNA = "SO:0000185"

PRODUCTS = {"CDS", "ncRNA", "tRNA", "rRNA"}

Table = ty.List[ty.Dict[str, str]]


def tables(handle) -> ty.Dict[str, Table]:
    found = {}
    for name, value in json.load(handle)["data"].items():
        lines = [l for l in value["content"].splitlines() if l and l[0] != "#"]
        # Headers are numbered: 1)id, 2)name, ...
        header = [c.split(")", 1)[1].strip() for c in lines[0].split("\t")]
        found[name] = [
            dict(zip(header, (v.strip() for v in l.split("\t")))) for l in lines[1:]
        ]
    return found


def units(tables) -> ty.Iterable[ty.Tuple]:
    promoters = {p["id"]: p for p in tables["promoters"]}
    terminators = {t["id"]: t for t in tables["terminators"]}
    for tu in tables["tus"]:
        promoter = promoters.get(tu.get("promoterId"))
        if tu.get("confidenceLevel") not in CONFIDENCE or not promoter:
            continue
        if not promoter.get("posTSS"):
            continue
        strand = 1 if promoter["strand"] == "forward" else -1
        tss = int(promoter["posTSS"])
        # A TU with several terminators makes one transcript for each
        for terminator_id in filter(None, tu.get("terminatorIds", "").split(",")):
            terminator = terminators[terminator_id]
            if strand == 1:
                start, stop = tss, int(terminator["rightEndPos"])
            else:
                start, stop = int(terminator["leftEndPos"]), tss
            if start >= stop:
                continue
            yield tu, f"{tu['id']}_{terminator_id}", start, stop, strand


def rna_type(inside) -> ty.Optional[str]:
    kinds = {f.type for f in inside}
    if kinds == {"tRNA"}:
        return PRECURSOR_RNA
    classes = {coding.first(f, "ncRNA_class") for f in inside}
    if kinds == {"ncRNA"} and len(classes) == 1:
        return as_so_term(classes.pop() or "ncRNA")
    return None


def parse(handle: ty.IO, genome: ty.IO) -> ty.Iterable[data.Entry]:
    found = tables(handle)
    record = next(SeqIO.parse(genome, "embl"))
    taxid = embl.taxid(record)
    named = collections.defaultdict(list)
    for feature in record.features:
        if feature.type in PRODUCTS and "pseudo" not in feature.qualifiers:
            for name in feature.qualifiers.get("gene", []):
                named[name].append(feature)

    for tu, accession, start, stop, strand in units(found):
        # RegulonDB names the genes; INSDC annotates some past RegulonDB's terminators
        genes = filter(None, (g.strip() for g in tu.get("tuGenes", "").split(";")))
        inside = [
            f
            for g in genes
            for f in named.get(g, [])
            if f.location.strand == strand
            and f.location.start < stop
            and f.location.end >= start
        ]
        pmids = dict.fromkeys(filter(None, tu.get("PMIDs", "").split(";")))
        references = [pubs.reference(p) for p in pmids]
        single = inside[0] if len(inside) == 1 else None
        transcript = coding.Transcript(
            accession=accession,
            spans=[(start, stop)],
            strand=strand,
            cds=None,
            name=f"{record.annotations['organism']} {tu['name']}",
            gene=coding.first(single, "gene") if single else None,
            locus_tag=coding.first(single, "locus_tag") if single else None,
            url=URL.format(operon=tu["operonId"]),
        )

        cds = [f for f in inside if f.type == "CDS"]
        if cds:
            # UTRs come from the genome: RegulonDB's UTR file is a base off on minus
            first = min(int(f.location.start) + 1 for f in cds)
            last = max(int(f.location.end) for f in cds)
            transcript = attr.evolve(transcript, cds=(first, last))
            yield from coding.entries(
                record, transcript, taxid, ASSEMBLY, "REGULONDB", references
            )
            continue

        so_term = rna_type(inside)
        if not so_term:
            continue
        label = SO_INSDC_MAPPING[so_term].replace("_", " ")
        entry = coding.as_entry(
            record,
            transcript,
            transcript.spans,
            accession,
            so_term,
            label,
            taxid,
            ASSEMBLY,
            "REGULONDB",
            references,
        )
        if entry.is_valid():
            yield entry
