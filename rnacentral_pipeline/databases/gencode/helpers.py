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

import attr

from rnacentral_pipeline.databases.helpers import publications as pub


def xref_data(entry):
    """
    Get GENCODE specific xref data. This will strip out the standard
    xref_data and add ENSEMBL as an xref.
    """

    xrefs = {k: v for k, v in entry.xref_data.items()}
    xrefs['Ensembl'] = [entry.accession]
    return xrefs


def primary_id(ensembl_entry):
    """
    Get a GENCODE primary id. This will use the ID prefixed with OTT and
    strip out the extra bits to create a GENCODE specific id.
    """
    return ensembl_entry.primary_id


def references():
    """
    Get the standard GENCODE reference.
    """
    return pub.reference(22955987)


def accession(entry):
    """
    This will produce an accession that is unique and specific to GENCODE.
    """
    return 'GENCODE:%s' % entry.accession


def update_entry(entry):
    """
    Modify an Ensembl Entry into a GENCODE Entry.
    """
    return attr.assoc(
        entry,
        accession=accession(entry),
        database='GENCODE',
        xref_data=xref_data(entry),
        optional_id=None,
        references=references(),
    )
