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

import typing as ty

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

LNC_DATABASES = {
    'lncipedia',
    'lncrnadb',
}


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

    if 'tRNA' in rna_type and sequence.length > 700:
        rna_type.remove('tRNA')
        return rna_type or set(['other'])

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


def prefer_lnc_over_anti(rna_type, _):
    """
    This will remove antisense_RNA to use lncRNA if we have both. This is
    because the term antisense_RNA usually but not always means
    antisense_lnc_RNA. Until we switch to SO terms which can express this we
    will switch.
    """

    if rna_type == set(['antisense_RNA', 'lncRNA']):
        return set(['lncRNA'])
    return rna_type


def remove_ncrna_if_possible(rna_types, _):
    """
    ncRNA is always consitent with everything else, so ignore it if possible.
    """

    if 'ncRNA' in rna_types and len(rna_types) > 1:
        rna_types.discard('ncRNA')
        return rna_types
    return rna_types


def from_lnc_database(data):
    """
    Check if this data came from any database that focuses on lncRNA's.
    """
    return any(a.database in LNC_DATABASES for a in data.accessions)


def from_ensembl(data):
    return any(a.database.lower() == 'ensembl' for a in data.accessions)


def choose_snoRNA_over_snRNA(rna_types, _):
    if rna_types == {'snoRNA', 'snRNA'}:
        return {'snoRNA'}
    return rna_types


def change_ncrna_to_other(rna_types, _):
    if rna_types == {'ncRNA'}:
        return {'other'}
    return rna_types


def remove_anti_if_required(rna_types, data):
    if 'antisense_RNA' in rna_types and from_lnc_database(data):
        rna_types.discard('antisense_RNA')
        rna_types.add('lncRNA')

    if 'antisense_RNA' in rna_types and from_ensembl(data):
        rna_types.discard('antisense_RNA')
        rna_types.add('lncRNA')

    return rna_types


def rna_type_of(data) -> ty.Optional[str]:
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
        raise ValueError("Could not find any database this sequence is from: %s", data)

    trusted = databases.intersection(TRUSTED_DATABASES)
    LOGGER.debug("Found %i trusted databases", len(trusted))
    if len(trusted) == 1:
        trusted = trusted.pop()
        rna_types = {acc.rna_type for acc in data.accessions}
        LOGGER.debug("Found %i rna_types from trusted dbs", len(rna_types))
        if len(rna_types) == 1:
            return rna_types.pop()

    rna_type = {acc.rna_type for acc in data.accessions}
    LOGGER.debug("Initial rna_types: %s", rna_type)

    corrections = [
        correct_other_vs_misc,
        remove_ambiguous,
        remove_ribozyme_if_possible,
        correct_by_length,
        prefer_lnc_over_anti,
        remove_ncrna_if_possible,
        change_ncrna_to_other,
        choose_snoRNA_over_snRNA,
        remove_anti_if_required,
    ]
    for correction in corrections:
        rna_type = correction(rna_type, data)
    rna_type.discard(None)

    LOGGER.debug("Corrected to %s rna_types", rna_type)
    if len(rna_type) == 1:
        LOGGER.debug("Found rna_type of %s for %s", rna_type, data.upi)
        return rna_type.pop()

    if not rna_type:
        LOGGER.debug("Corrections removed all rna_types")
    LOGGER.debug("Using fallback count method for %s", data.upi)

    counts = Counter(a.rna_type for a in data.accessions if a.rna_type in rna_type)
    selected = counts.most_common(1)[0][0]
    if not selected:
        raise ValueError("Fallback selection failed")
    return selected
