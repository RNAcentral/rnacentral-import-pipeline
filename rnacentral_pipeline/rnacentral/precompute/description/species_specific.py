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

import re
import logging
import operator as op
import itertools as it

from .. import utils

LOGGER = logging.getLogger(__name__)


CHOICES = {
    'miRNA': ['mirbase', 'refseq', 'hgnc', 'gencode', 'rfam', 'ensembl', 'ena'],
    'precursor_RNA': ['mirbase', 'refseq', 'rfam', 'hgnc', 'gencode', 'ensembl', 'ena'],
    'ribozyme': ['refseq', 'rfam', 'pdbe', 'ensembl', 'ena'],
    'hammerhead_ribozyme': ['refseq', 'rfam', 'pdbe', 'ensembl', 'ena'],
    'autocatalytically_spliced_intron': ['refseq', 'rfam', 'pdbe', 'ensembl', 'ena'],

    '__generic__': [
        'mirbase',
        'wormbase',
        'hgnc',
        'gencode',
        'ensembl',
        'tair',
        'sgd',
        'flybase',
        'dictbase',
        'pombase',
        'mgi',
        'rgd',
        'lncipedia',
        'lncrnadb',
        'gtrnadb',
        'tmrna website',
        'pdbe',
        'refseq',
        'rfam',
        'modomics',
        'vega',
        'srpdb',
        'snopy',
        'silva',
        'greengenes',
        'rdp',
        'ena',
        'noncode',
    ]
}

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
    if required_rna_type in set(['miRNA', 'precursor_RNA']):
        allowed_rna_types = set(['miRNA', 'precursor_RNA'])
    allowed_rna_types.add('ncRNA')

    def fn(db_name, accession):
        if accession.database != db_name:
            return False

        # PDBe has lots of things called 'misc_RNA' that have a good
        # description, so we allow this to use PDBe's misc_RNA descriptions
        if accession.database == 'pdbe' and accession.rna_type == 'misc_rna':
            return True

        return accession.rna_type in allowed_rna_types
    return fn


def accept_any(db_name, accession):
    return accession.database == db_name


def description_order(name):
    """
    Computes a tuple to order descriptions by.
    """
    return (round(utils.entropy(name), 3), [-ord(c) for c in name])


def select_best_description(descriptions):
    """
    This will generically select the best description. We select the string
    with the maximum entropy and lowest description. The entropy constraint is
    meant to deal with names for PDBe which include things like AP*CP*... and
    other repetitive databases. The other constraint is to try to select things
    that come from a lower number (if numbered) item.
    """
    return max(descriptions, key=description_order)


def compute_item_ranges(items):
    names = []
    data = sorted(utils.item_sorter(item) for item in items if item)
    for item, numbers in it.groupby(data, op.itemgetter(0)):
        if not item:
            continue

        range_format = '%i-%i'
        if '-' in item:
            range_format = ' %i to %i'
        for (start, stop) in utils.group_consecutives(n[1] for n in numbers):
            if stop is None:
                names.append(item + str(start))
            else:
                prefix = item
                if prefix.endswith('-'):
                    prefix = prefix[:-1]
                names.append(prefix + range_format % (start, stop))

    return names


def add_term_suffix(base, additional_terms, name, max_items=3):
    items = compute_item_ranges(additional_terms)

    suffix = 'multiple %s' % name
    if len(items) < max_items:
        suffix = ', '.join(items)

    if suffix in base:
        return base

    return '{basic} ({suffix})'.format(
        basic=base.strip(),
        suffix=suffix,
    )


def select_with_several_genes(accessions, name, pattern,
                              description_items=None,
                              attribute='gene',
                              max_items=3):
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
                description += ' (%s)' % suffix
        return description

    regexp = pattern % getter(candidate)
    basic = re.sub(regexp, '', candidate.description)

    func = getter
    if description_items is not None:
        func = op.attrgetter(description_items)

    items = sorted([func(a) for a in accessions if func(a)], key=utils.item_sorter)
    if not items:
        return basic

    return add_term_suffix(basic, items, name, max_items=max_items)


def improve_predicted_description(rna_type, accessions, description):
    alt = []
    for accession in accessions:
        if accession.database == 'rfam' and \
                accession.rna_type == rna_type:
            alt.append(accession)

    if not alt:
        return description

    description = select_best_description([a.description for a in alt])

    # If there is a gene to append we should
    genes = [acc.gene for acc in accessions]
    if genes:
        description = add_term_suffix(description, genes, 'genes')

    return description


def cleanup(rna_type, db_name, description):
    # There are often some extra terms we need to strip
    description = utils.remove_extra_description_terms(description)
    if db_name == 'RefSeq':
        description = utils.trim_trailing_rna_type(rna_type, description)

    return description.strip()


def improve_mirbase_description(rna_type, accessions):
    product_name = 'precursors'
    if rna_type == 'miRNA':
        product_name = 'miRNAs'

    return select_with_several_genes(
        accessions,
        product_name,
        r'\w+-%s',
        description_items='optional_id',
        max_items=5)


def is_mod(db_name):
    """
    Check if this database is a MOD that we can generate mulitple gene
    descriptions.
    """
    return db_name in set(['HGNC', 'Ensembl', 'GENCODE'])


def improve_mod_description(accessions):
    """
    Improve the description from a model organism database (MOD).
    """
    return select_with_several_genes(
        accessions,
        'genes',
        r'\(%s\)$',
        max_items=5)


def replace_nulls(rna_type, description):
    if rna_type not in description:
        return re.sub('null$', rna_type, description)
    return description


def description_of(rna_type, sequence):
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

    if rna_type not in CHOICES:
        LOGGER.debug("Falling back to generic ordering for %s", rna_type)

    ordering = CHOICES.get(rna_type, CHOICES['__generic__'])
    selector = suitable_xref(rna_type)

    try:
        db_name, accessions = utils.best(ordering, sequence.accessions, selector)
    except utils.NoBestFoundException:
        db_name, accessions = utils.best(ordering, sequence.accessions, accept_any)

    # It is possible that one sequence maps to several HGNC genes and we should
    # indicate this
    if is_mod(db_name):
        description = improve_mod_description(accessions)

    # Similar issue for miRBase sequences
    elif db_name == 'miRBase':
        description = improve_mirbase_description(rna_type, accessions)

    # Fall back to a simple generic method
    else:
        descriptions = [accession.description for accession in accessions]
        description = select_best_description(descriptions)

    # Sometimes we get a description that is 'predicted' from some databases.
    # It would be better to pull from Rfam which may have a more useful
    # description.
    if 'predicted' in description:
        description = improve_predicted_description(
            rna_type,
            accessions,
            description,
        )

    # ENA sometimes has things that end with 'null', which is bad.
    if description.endswith(' null'):
        description = replace_nulls(rna_type, description)

    return cleanup(rna_type, db_name, description)
