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
from collections import Counter

LOGGER = logging.getLogger(__name__)

TRUSTED_DATABASES = set([
    'miRBase'
])
"""This defines the set of databases that we always trust. If we have
annotations from one of these databases we will always use if, assuming all
annotations from the database agree.
"""


def correct_by_length(rna_type, sequence):
    """
    This will correct the miRNA/precursor_RNA conflict and ambiguitity. Some
    databases like 'HGNC' will call a precursor_RNA miRNA. We correct this
    using the length of the sequence as well as using the length to distinguish
    between the two.
    """

    if rna_type == set([u'precursor_RNA', u'miRNA']) or \
            rna_type == set(['miRNA']):
        if 15 <= sequence.length <= 30:
            return set(['miRNA'])
        return set(['precursor_RNA'])
    return rna_type


def correct_other_vs_misc(rna_type, _):
    """
    Given 'misc_RNA' and 'other' we prefer 'other' as it is more specific. This
    will only select 'other' if 'misc_RNA' and other are the only two current
    rna_types.
    """

    if rna_type == set(['other', 'misc_RNA']):
        return set(['other'])
    return rna_type


def remove_ambiguous(rna_type, _):
    """
    If there is an annotation that is more specific than other or misc_RNA we
    should use it. This will remove the annotations if possible
    """

    ambiguous = set(['other', 'misc_RNA'])
    specific = rna_type.difference(ambiguous)
    if specific:
        return specific
    return rna_type


def remove_ribozyme_if_possible(rna_type, _):
    """
    This will remove the ribozyme rna_type from the set of rna_types if there
    is a more specific ribozyme annotation avaiable.
    """

    ribozymes = set(['hammerhead', 'hammerhead_ribozyme',
                     'autocatalytically_spliced_intron'])
    if 'ribozyme' in rna_type and rna_type.intersection(ribozymes):
        rna_type.discard('ribozyme')
        return rna_type
    return rna_type


def remove_ncrna_if_possible(rna_types, _):
    """
    ncRNA is always consitent with everything else, so ignore it if possible.
    """

    if 'ncRNA' in rna_types and len(rna_types) > 1:
        rna_types.discard('ncRNA')
        return rna_types
    return rna_types


def get_rna_types_from(xrefs, name):
    """
    Determine the rna_types as annotated by some database.

    Parameters
    ----------
    xrefs : iterable
        The list of xrefs to fitler to extract the rna types from.
    name : str
        The name of the database to use.

    Returns
    -------
    rna_types : set
        A set of rna types that is annotated by the given database.
    """

    rna_types = set()
    for xref in xrefs:
        if xref.db.name == name:
            rna_types.add(xref.accession.get_rna_type())
    return rna_types


def rna_type_of(data):
    """
    Determine the rna_type for a given sequence and collection of xrefs. The
    idea behind this is that not all databases are equally trustworthy, so some
    annotations of rna_type should be ignored while others should be trusted.
    The goal of this function then is to examine all annotated rna_types and
    selected the annotation(s) that are most reliable for the given sequence.

    Note that if this cannot use any of the normal rules to select a single
    rna_type then it will use a simple method of taking the most common,
    followed by alphabetically first rna_type.

    Parameters
    ----------
    sequence : Rna
        The sequence to examine.
    xrefs : iterable
        The collection of xrefs to use.

    Returns
    -------
    rna_type : str
        The selected RNA type for this sequence.
    """

    databases = {acc.database for acc in data.accessions}
    if not databases:
        LOGGER.error("Could not find any database this sequence is from: %s",
                     data)
        return None

    trusted = databases.intersection(TRUSTED_DATABASES)
    LOGGER.debug("Found %i trusted databases", len(trusted))
    if len(trusted) == 1:
        trusted = trusted.pop()
        rna_types = get_rna_types_from(xrefs, trusted)
        LOGGER.debug("Found %i rna_types from trusted dbs", len(rna_types))
        if len(rna_types) == 1:
            return rna_types.pop()

    accessions = []
    rna_type = set()
    for accession in data.accessions:
        rna_type.add(accession.rna_type)
        accessions.append(accession)

    LOGGER.debug("Initial rna_types: %s", rna_type)

    corrections = [
        correct_other_vs_misc,
        remove_ambiguous,
        remove_ribozyme_if_possible,
        remove_ncrna_if_possible,
        correct_by_length,
    ]
    for correction in corrections:
        rna_type = correction(rna_type, sequence)

    LOGGER.debug("Corrected to %s rna_types", rna_type)
    if len(rna_type) == 1:
        LOGGER.debug("Found rna_type of %s for %s", rna_type, sequence.upi)
        return rna_type.pop()

    if not rna_type:
        LOGGER.debug("Corrections removed all rna_types")
    LOGGER.debug("Using fallback count method for %s", sequence.upi)

    counts = Counter(x.accession.get_rna_type() for x in xrefs)
    return counts.most_common(1)[0][0]
