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
import math
import string
import logging
import operator as op
import itertools as it
from collections import Counter

LOGGER = logging.getLogger(__name__)


class NoBestFoundException(Exception):
    """
    Raised when it is not possible to find the 'best' xref.
    """

    pass


def best(ordered_choices, possible, check, default=None):
    """
    Select the best xref from several possible xrefs given the ordered list of
    xref database names. This function will iterate over each database name and
    select all xrefs that come from the first (most preferred) database. This
    uses the check function to see if the database contains the correct
    information because this doesn't always just check based upon the database
    names or xref, but also the rna_type (in some cases). Using a function
    gives a lot of flexibility in how we select the acceptable xrefs.

    Parameters
    ----------
    ordered_choices : list
        A list of several possible xref database names. These should be in the
        order in which they are preferred.
    possible : list
        The list of xrefs to find the best for.
    check : callable
        A callable object to see if given xref and database name match.
    default : obj, None
        The default value to return if we cannot find a good xref.

    Returns
    -------
    selected : obj
        The list of xrefs which are 'best' given the choices. If there is no
        good xref the default value is returned.
    """

    for choice in ordered_choices:
        found = [entry for entry in possible if check(choice, entry)]
        if found:
            return (choice, found)
    raise NoBestFoundException("Could not select a best for: %s" % str(possible))


def entropy(data):
    """
    This computes an approximation of the entropy in the given string. This is
    used to select the 'best' name of several possible names. The reason we use
    this and not just length is that in some cases (like sequences from
    structures) the name will be very long because it contains the sequence
    itself. For example:

    RNA (5'-R(*GP*UP*GP*GP*UP*CP*UP*GP*AP*UP*GP*AP*GP*GP*CP*C)-3') from synthetic construct (PDB 3D0M, chain X)

    This is not a useful name, but it is very long. Thus we do not want it.
    What we are generally after is something with the most information (to a
    human) in it. We get at this by using entropy as a guide. To make sure this
    computes a useful entropy we restrict the string to only the letters and
    numbers while ignoring case. This should provide a clean measure of an
    informative string for a person.

    Parameters
    ----------
    data : str
        A string to get the entropy for.

    Returns
    -------
    entropy : float
        The entropy of the string.
    """

    if len(data) <= 1:
        return 0
    allowed_characters = set(string.ascii_lowercase + string.digits)
    valid = list(d.lower() for d in data if d.lower() in allowed_characters)
    counts = Counter(valid)
    probs = [float(c) / len(valid) for c in counts.values() if c > 0]
    return sum(-1 * p * math.log(p, 2) for p in probs)


def item_sorter(name):
    match = re.search(r"(\d+)$", name)
    if match:
        name = re.sub(r"(\d+)$", "", name)
        return (name, int(match.group(1)))
    return (name, "")


def group_consecutives(data, min_size=2):
    """
    Modified from the python itertools docs.
    """

    def grouper(ix):
        if isinstance(ix[0], int) and isinstance(ix[1], int):
            return ix[0] - ix[1]
        else:
            return float("inf")

    for _, group in it.groupby(enumerate(data), grouper):
        key = op.itemgetter(1)
        group = [key(g) for g in group]
        if len(group) > 1:
            if group[-1] - group[0] < min_size:
                for member in group:
                    yield member, None
            else:
                yield group[0], group[-1]
        else:
            yield group[0], None


def remove_extra_description_terms(description):
    """
    Sometimes the description contains some things we don't want to include
    like trailing '.', extra spaces, the string (non-protein coding) and a typo
    in the name of tmRNA's. This corrects those issues.
    """

    description = re.sub(r"\(\s*non\s*-\s*protein\s+coding\s*\)", "", description)
    description = re.sub(
        r"transfer-messenger mRNA", "transfer-messenger RNA", description, re.IGNORECASE
    )
    description = re.sub(r"\s\s+", " ", description)
    description = re.sub(r"\.?\s*$", "", description)
    return description


def trim_trailing_rna_type(rna_type, description):
    """
    Some descriptions like the ones from RefSeq (URS0000ABD871/9606) and some
    from flybase (URS0000002CB3/7227) end with their RNA type. This isn't
    really all that useful as we display the RNA type anyway. So  this
    out along with ', ' that tends to go along with those.
    """

    if "predicted gene" in description:
        return description

    trailing = [rna_type]
    if rna_type == "lncRNA" or "antisense" in rna_type:
        trailing.append("long non-coding RNA")

    if rna_type == "telomerase_RNA":
        trailing.append("telomerase RNA")

    for value in trailing:
        pattern = r"[, ]*%s$" % value
        description = re.sub(pattern, "", description, re.IGNORECASE)
    return description
