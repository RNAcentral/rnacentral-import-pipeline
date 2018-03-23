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


import attr

URL = 'https://www.ncbi.nlm.nih.gov/nuccore/{primary_id}.{version}'


def optional_id(entry):
    key = 'GeneID'
    gene_ids = entry.xref_data[key]
    if not gene_ids:
        return None
    return key + ':%s' % gene_ids[0]


def external_id(entry):
    """
    Find the RefSeq external id.
    """

    ena_refs = entry.xref_data['ena_refs']
    if not ena_refs:
        raise ValueError("Must have ena_refs for %s" % entry)
    result, _ = ena_refs['REFSEQ']
    return result


def general_accession(entry):
    return entry.accession.split(':')[0].split('.')


def seq_version(entry):
    return general_accession(entry)[-1]


def parent_accession(entry):
    return general_accession(entry)[0]


def url(entry):
    return URL.format(
        primary_id=external_id(entry),
        version=seq_version(entry),
    )


def as_entry(entry):
    """
    Modify an ENA entry into an approbate RefSeq entry.
    """
    return attr.evolve(
        entry,
        database='REFSEQ',
        exons=[],
        seq_version=seq_version(entry),
        parent_accession=parent_accession(entry),
        url=url(entry),
        primary_id=external_id(entry),
        optional_id=optional_id(entry),
    )
