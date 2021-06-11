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

import itertools as it
import logging
import operator as op
import re
import typing as ty

import attr

from rnacentral_pipeline.databases.data import Database, RnaType
from rnacentral_pipeline.databases.sequence_ontology import tree
from rnacentral_pipeline.rnacentral.precompute.data import context
from rnacentral_pipeline.rnacentral.precompute.data import sequence as seq
from rnacentral_pipeline.rnacentral.precompute.data.accession import Accession
from rnacentral_pipeline.rnacentral.precompute.qa import contamination as cont

from rnacentral_pipeline.rnacentral.precompute import utils

LOGGER = logging.getLogger(__name__)


ORDERING = [
    Database.mirbase,
    Database.wormbase,
    Database.hgnc,
    Database.gencode,
    Database.ensembl,
    Database.tair,
    Database.sgd,
    Database.flybase,
    Database.dictybase,
    Database.pombase,
    Database.mgi,
    Database.rgd,
    Database.zfin,
    Database.mirgenedb,
    Database.lncipedia,
    Database.lncrnadb,
    Database.lncbook,
    Database.gtrnadb,
    Database.tmrna_website,
    Database.five_srrnadb,
    Database.pdbe,
    Database.refseq,
    Database.ensembl_plants,
    Database.ensembl_metazoa,
    Database.ensembl_protists,
    Database.ensembl_fungi,
    Database.genecards,
    Database.malacards,
    Database.intact,
    Database.rfam,
    Database.tarbase,
    Database.lncbase,
    Database.snodb,
    Database.snorna_database,
    Database.pirbase,
    Database.modomics,
    Database.vega,
    Database.srpdb,
    Database.snopy,
    Database.crw,
    Database.silva,
    Database.greengenes,
    Database.rdp,
    Database.ena,
    Database.zwd,
    Database.noncode,
]
"""
A dict that defines the ordered choices for each type of RNA. This is the
basis of our name selection for the rule based approach. The fallback,
__generic__ is a list of all database roughly ordered by how good the names
from each one are.
"""


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


class DatabaseSpecifcNameBuilder(object):
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


def suitable_xref(required_rna_type):
    """
    Create a function, based upon the given rna_type, which will test if the
    database has assinged the  correct rna_type to the sequence. This is used
    when selecting the description so we use the description from a database
    that gets the rna_type correct. There are exceptions for
    miRNA/precursor_RNA as well as PDBe's misc_RNA information.

    Parameters
    ----------
    required_rna_type : str
        The rna_type to use to check the database.

    Returns
    -------
    fn : function
        A function to detect if the given xref has information about the
        rna_type that can be used for determining the rna_type and description.
    """

    # allowed_rna_types is the set of rna_types which a database is allowed to
    # call the sequence for this function to trust the database's opinion on
    # the description/rna_type. We allow database to use one of
    # miRNA/precursor_RNA since some databases (Rfam, HGNC) do not correctly
    # distinguish the two but do have good descriptions otherwise.
    allowed_rna_types = set([required_rna_type])
    if required_rna_type in set(["miRNA", "precursor_RNA"]):
        allowed_rna_types = set(["miRNA", "precursor_RNA"])
    allowed_rna_types.add("ncRNA")

    def fn(db_name, accession):
        if accession.database != db_name:
            return False

        # PDBe has lots of things called 'misc_RNA' that have a good
        # description, so we allow this to use PDBe's misc_RNA descriptions
        if accession.database is Database.pdbe and accession.rna_type == "misc_rna":
            return True

        return accession.rna_type in allowed_rna_types

    return fn


def accept_any(db_name: str, accession: Accession) -> bool:
    return accession.database == db_name


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

        for (start, stop) in utils.group_consecutives(n[1] for n in numbers):
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
    candidate = min(accessions, key=getter)
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


def improve_predicted_description(
    rna_type: str, accessions: ty.List[Accession], description: str
) -> str:
    alt = []
    for accession in accessions:
        if accession.database == "rfam" and accession.rna_type == rna_type:
            alt.append(accession)

    if not alt:
        return description

    description = select_best_description([a.description for a in alt])

    # If there is a gene to append we should
    genes = [acc.gene for acc in accessions]
    if genes:
        description = add_term_suffix(description, genes, "genes")

    return description


def cleanup(rna_type: str, db_name: str, description: str) -> str:
    # There are often some extra terms we need to strip
    description = utils.remove_extra_description_terms(description)
    description = description.replace("()", "")
    if db_name == "refseq":
        description = utils.trim_trailing_rna_type(rna_type, description)

    if db_name == "tarbase":
        description = description.replace("TARBASE:", "")
    description = description.replace(" (None)", "")
    description = re.sub(r"\s\s+", " ", description)
    return description.strip()


def replace_nulls(rna_type: str, description: str) -> str:
    if rna_type not in description:
        return re.sub("null$", rna_type, description)
    return description


def description_of(rna_type: str, sequence: seq.Sequence) -> str:
    """
    Determine the name for the species specific sequence. This will examine
    all descriptions in the xrefs and select one that is the 'best' name for
    the molecule. If no xref can be selected as a name, then None is returned.
    This can occur when no xref in the given iterable has a matching rna_type.
    The best description will be the one from the xref which agrees with the
    computed rna_type and has the maximum entropy as estimated by `entropy`.
    The reason this is used over length is that some descriptions which come
    from PDBe import are highly repetitive because they are for short sequences
    and they contain the sequence in the name. Using entropy basis away from
    those sequences to things that are hopefully more informative.

    Parameters
    ----------
    rna_type : str
        The type for the sequence

    sequence : Rna
        The sequence entry we are trying to select a name for

    xrefs : iterable
        An iterable of the Xref entries that are specific to a species for this
        sequence.

    Returns
    -------
    name : str, None
        A string that is a description of the sequence, or None if no sequence
        could be selected.
    """

    selector = suitable_xref(rna_type)
    try:
        db_name, accessions = utils.best(ORDERING, sequence.accessions, selector)
    except utils.NoBestFoundException:
        db_name, accessions = utils.best(ORDERING, sequence.accessions, accept_any)

    builder = DatabaseSpecifcNameBuilder()
    description = builder(db_name, rna_type, accessions)

    # Sometimes we get a description that is 'predicted' from some databases.
    # It would be better to pull from Rfam which may have a more useful
    # description.
    if "predicted" in description:
        description = improve_predicted_description(
            rna_type,
            accessions,
            description,
        )

    # ENA sometimes has things that end with 'null', which is bad.
    if description.endswith(" null"):
        description = replace_nulls(rna_type, description)

    result = cleanup(rna_type, db_name, description)
    if isinstance(result, bytes):
        return result.decode("utf-8")
    elif isinstance(result, str):
        return result
    else:
        raise ValueError("Unknown type of description")
