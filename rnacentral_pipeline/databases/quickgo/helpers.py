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
import operator as op

import databases.helpers.publications as pub

from . import data

go_id = op.itemgetter('GO_ID')
evidence = op.itemgetter('Evidence')
assigned_by = op.itemgetter('Assigned_by')


def rna_id(entry):
    """
    Get the UPI for the entry, or fail if there is none.
    """

    if entry['DB'] == 'RNAcentral':
        return entry['DB_Object_ID']
    raise ValueError("All entries are expected to come from RNAcentral")


def qualifier(entry):
    """
    Get the qualifer for this entry.
    """
    return entry['Qualifier'][0]


def publications(entry):
    references = []
    for reference in entry['DB:Reference']:
        match = re.match(r'^PMID:(\d+)$', reference)
        if match:
            references.append(pub.reference(match.group(1)))
    return references


def extensions(record):
    result = []
    for extension in record['Annotation Extension']:
        for part in extension.split(','):
            match = re.match(r'(\w+)\((.+)\)', part)
            if match:
                result.append(data.AnnotationExtension(
                    qualifier=match.group(1),
                    target=match.group(2),
                ))
    return result
