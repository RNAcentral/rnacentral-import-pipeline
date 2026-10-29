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

import logging
import re
import typing as ty

from rnacentral_pipeline.databases.data import Database, RnaType
from rnacentral_pipeline.databases.sequence_ontology import tree
from rnacentral_pipeline.rnacentral.precompute import utils
from rnacentral_pipeline.rnacentral.precompute.data import context
from rnacentral_pipeline.rnacentral.precompute.data import sequence as seq
from rnacentral_pipeline.rnacentral.precompute.data.accession import Accession
from rnacentral_pipeline.rnacentral.precompute.description.name_builders import (
    DatabaseSpecificNameBuilder,
    add_term_suffix,
    select_best_description,
)
from rnacentral_pipeline.rnacentral.precompute.qa import contamination as cont

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
    Database.plncdb,
    Database.lncipedia,
    Database.lncrnadb,
    Database.lncbook,
    Database.gtrnadb,
    Database.tmrna_website,
    Database.five_srrnadb,
    Database.ribocentre,
    Database.pdbe,
    Database.refseq,
    Database.ensembl_plants,
    Database.ensembl_metazoa,
    Database.ensembl_protists,
    Database.ensembl_fungi,
    Database.genecards,
    Database.malacards,
    Database.intact,
    Database.expression_atlas,
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
    Database.evlncrnas,
    Database.psicquic,
    Database.ribovision,
    Database.mirtrondb,
    Database.japonicusdb,
    Database.circatlas,
    Database.circpedia,
    Database.mgnify,
]
"""
A dict that defines the ordered choices for each type of RNA. This is the
basis of our name selection for the rule based approach. The fallback,
__generic__ is a list of all database roughly ordered by how good the names
from each one are.
"""


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
    the molecule. Every RNA must get a description, so if no xref at all can
    be selected as a name (e.g. a since-superseded sequence with no active
    accessions), a NoBestFoundException is raised instead of returning None.
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
    name : str
        A string that is a description of the sequence.
    """

    selector = suitable_xref(rna_type)
    try:
        db_name, accessions = utils.best(ORDERING, sequence.accessions, selector)
    except utils.NoBestFoundException:
        # Every RNA must get a description, so if there is no active
        # accession at all (e.g. a since-superseded sequence) to fall back
        # to, this propagates and fails the pipeline rather than silently
        # producing no description.
        db_name, accessions = utils.best(ORDERING, sequence.accessions, accept_any)

    builder = DatabaseSpecificNameBuilder()
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

    return cleanup(rna_type, db_name, description)
