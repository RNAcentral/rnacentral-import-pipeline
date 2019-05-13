# -*- coding: utf-8 -*-

"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
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


def short_description(description, sequence):
    """
    Shorten a description to remove the leading species (common name) from
    description. This is useful for display in the genome browser.
    """

    common = set()
    species = set()
    for acc in sequence.accessions:
        species.update(acc.all_species)
        common.update(u'(%s)' % c for c in acc.all_common_names)
    leading = [l for l in list(species) + list(common) if l]

    for name in leading:
        # This is used instead of building a regex to avoid dealing with issues
        # due to parenthesis.
        name = name
        if description.startswith(name):
            description = description[len(name):]
        description = description.strip()

    description = re.sub(r'^\s*(\(\))?\s*', '', description)

    if description.startswith('(') and description.endswith(')'):
        description = description[1:-1]

    return description
