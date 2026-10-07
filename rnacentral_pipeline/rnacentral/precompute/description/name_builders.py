# -*- coding: utf-8 -*-

"""
Copyright [2009-2017] EMBL-European Bioinformatics Institute
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

"""
Per-database description building: given a set of accessions that all agree
on which database and rna_type to use (selected by
species_specific.description_of), work out the best human-readable
description for them. Split out of species_specific.py since it's a
self-contained concern (name selection) distinct from that module's job of
choosing *which* database's accessions to use in the first place.
"""

import itertools as it
import operator as op
import re
import typing as ty

import attr

from rnacentral_pipeline.databases.data import Database
from rnacentral_pipeline.rnacentral.precompute import utils
from rnacentral_pipeline.rnacentral.precompute.data.accession import Accession


def description_order(name: str):
    """
    Computes a tuple to order descriptions by.
    """
    return (round(utils.entropy(name), 3), [-ord(c) for c in name])


def select_best_description(descriptions: ty.List[str]):
    """
    This will generically select the best description. We select the string
    with the maximum entropy and lowest description. The entropy constraint is
    meant to deal with names for PDBe which include things like AP*CP*... and
    other repetitive databases. The other constraint is to try to select things
    that come from a lower number (if numbered) item.
    """
    return max(descriptions, key=description_order)


def compute_item_ranges(items):
    data = sorted(utils.item_sorter(item) for item in items if item)
    grouped = it.groupby(data, op.itemgetter(0))
    names = []
    for gene, numbers in grouped:
        if not gene:
            continue

        range_format = "%i-%i"
        if "-" in gene:
            range_format = " %i to %i"

        for start, stop in utils.group_consecutives(n[1] for n in numbers):
            if stop is None:
                names.append(gene + str(start))
            else:
                prefix = gene
                if prefix.endswith("-"):
                    prefix = prefix[:-1]
                names.append(prefix + range_format % (start, stop))

    return names


def add_term_suffix(base, additional_terms, name: str, max_items=3):
    items = compute_item_ranges(additional_terms)

    suffix = "multiple %s" % name
    if len(items) < max_items:
        suffix = ", ".join(items)

    if suffix in base:
        return base

    return "{basic} ({suffix})".format(
        basic=base.strip(),
        suffix=suffix,
    )


def select_with_several_genes(
    accessions: ty.List[Accession],
    name: str,
    pattern,
    description_items=None,
    attribute="gene",
    max_items=3,
):
    """
    This will select the best description for databases where more than one
    gene (or other attribute) map to a single URS. The idea is that if there
    are several genes we should use the lowest one (RNA5S1, over RNA5S17) and
    show the names of genes, if possible. This will list the genes if there are
    few, otherwise provide a note that there are several.
    """

    getter = op.attrgetter(attribute)
    # Sort accessions missing the attribute last: several can share a null
    # gene/locus_tag, and min() cannot compare None with None. Doing so also
    # keeps the pattern below off a candidate whose value is None.
    candidate = min(accessions, key=lambda a: (getter(a) is None, getter(a) or ""))
    genes = set(getter(a) for a in accessions if getter(a))
    if not genes or len(genes) == 1:
        description = candidate.description
        # Append gene name if it exists and is not present in the description
        # already
        if genes:
            suffix = genes.pop()
            if suffix not in description:
                description += " (%s)" % suffix
        return description

    regexp = pattern % getter(candidate)
    basic = re.sub(regexp, "", candidate.description)

    func = getter
    if description_items is not None:
        func = op.attrgetter(description_items)

    possible = {func(a) for a in accessions if func(a)}
    items = sorted(possible, key=utils.item_sorter)
    if not items:
        return basic

    return add_term_suffix(basic, items, name, max_items=max_items)


class DatabaseSpecificNameBuilder(object):
    def ensembl(self, accessions: ty.List[Accession], _) -> str:
        return select_with_several_genes(
            accessions, "genes", r"\(%s\)$", description_items="gene", max_items=5
        )

    def flybase(self, accessions: ty.List[Accession], _) -> str:
        return select_with_several_genes(
            accessions,
            "genes",
            r"%s$",
            attribute="locus_tag",
            description_items="locus_tag",
            max_items=6,
        )

    def gencode(self, accessions: ty.List[Accession], _) -> str:
        return select_with_several_genes(
            accessions, "genes", r"\(%s\)$", description_items="gene", max_items=5
        )

    def gtrnadb(self, accessions: ty.List[Accession], _) -> str:
        return select_with_several_genes(
            accessions, "tRNAs", r"\(%s\)$", description_items="gene", max_items=5
        )

    def hgnc(self, accessions: ty.List[Accession], _) -> str:
        return select_with_several_genes(
            accessions, "genes", r"\(%s\)$", description_items="gene", max_items=5
        )

    def mgi(self, accessions: ty.List[Accession], _) -> str:
        return select_with_several_genes(
            accessions, "genes", r"\(%s\)$", description_items="gene", max_items=5
        )

    def mirbase(self, accessions: ty.List[Accession], rna_type: str) -> str:
        if rna_type == "miRNA":
            return select_with_several_genes(
                accessions,
                "miRNAs",
                r"\w+-%s",
                description_items="optional_id",
                max_items=5,
            )

        updated = []
        for accession in accessions:
            gene = accession.optional_id
            if not gene and accession.description.endswith("stem-loop"):
                gene = accession.description.split(" ")[-2]
            if not gene and accession.description.endswith("stem loop"):
                gene = accession.description.split(" ")[-3]
            if not gene:
                last = accession.description.split(" ")[-1]
                if last.endswith("-3p") or last.endswith("-5p"):
                    last = last[:-4]
                if re.match(r"^.*-mir-\d+$", last, re.IGNORECASE):
                    gene = last
            if not gene:
                raise ValueError(f"Could not find gene for mirbase {accession}")
            match = re.match(r"^([^-]+?-mir-[^-]+)(.+)?$", gene)
            name = accession.description
            if match:
                full = match.group(1)
                parts = full.split("-", 3)
                trimmed = "-".join(parts[:3])
                name = "{species} ({common_name}) microRNA {gene} precursor".format(
                    species=accession.species,
                    common_name=accession.common_name or "",
                    gene=trimmed,
                )
            changed = attr.evolve(
                accession,
                description=name,
                gene=gene,
            )
            updated.append(changed)

        return select_with_several_genes(
            updated,
            "precursors",
            r"\w+-%s",
            description_items="optional_id",
            max_items=5,
        )

    def sgd(self, accessions: ty.List[Accession], _):
        return select_with_several_genes(
            accessions,
            "genes",
            r"%s$",
            attribute="optional_id",
            description_items="optional_id",
            max_items=6,
        )

    def tair(self, accessions: ty.List[Accession], _):
        return select_with_several_genes(
            accessions,
            "genes",
            r"%s$",
            attribute="locus_tag",
            description_items="locus_tag",
            max_items=6,
        )

    def _fallback(self, accessions: ty.List[Accession], _):
        descriptions = [accession.description for accession in accessions]
        return select_best_description(descriptions)

    def __call__(
        self, database: Database, rna_type: str, accessions: ty.List[Accession]
    ):
        name = database.normalized().lower()
        method = getattr(self, name, self._fallback)
        return method(accessions, rna_type)
