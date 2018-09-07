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

    species = {acc.species for acc in sequence.accessions if acc.species}
    common = {'(%s)' % acc.common_name for acc in sequence.accessions if acc.common_name}
    leading = list(species) + list(common)

    for name in leading:
        # This is used instead of building a regex to avoid dealing with issues
        # due to parenthesis.
        if description.startswith(name):
            description = description[len(name):]
        description = description.strip()

    if description.startswith('(') and description.endswith(')'):
        description = description[1:-1]
    description = re.sub(r'^\s+', '', description)
    return description
