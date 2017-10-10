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

from datetime import date
import operator as op
import itertools as it

import xml.etree.cElementTree as ET

from .data import XmlEntry


BASE_SQL = """
SELECT
    rna.upi,
    rna.md5,
    xref.taxid,
    xref.deleted,
    acc.species,
    acc.organelle organelles,
    acc.external_id,
    acc.optional_id,
    acc.non_coding_id,
    acc.accession,
    acc.function,
    acc.gene,
    acc.gene_synonym,
    acc.feature_name,
    acc.ncrna_class,
    acc.product,
    acc.common_name,
    acc.note,
    acc.parent_ac || '.' || acc.seq_version as parent_accession,
    db.display_name as expert_db,
    release1.timestamp as created,
    release2.timestamp as last,
    rna.len as length,
    pre.rna_type,
    pre.description,
    acc.locus_tag,
    acc.standard_name,
    ref.id pub_id,
    ref.location,
    ref.title pub_title,
    ref.pmid pubmed,
    ref.doi,
    ref.authors author,
    acc.classification tax_string
FROM
    xref,
    rnc_accessions acc,
    rnc_database db,
    rnc_release release1,
    rnc_release release2,
    rna,
    rnc_rna_precomputed pre,
    rnc_references ref,
    rnc_reference_map ref_map
WHERE
    xref.ac = acc.accession AND
    xref.dbid = db.id AND
    xref.created = release1.id AND
    xref.last = release2.id AND
    xref.upi = rna.upi AND
    xref.upi = pre.upi AND
    xref.taxid = pre.taxid AND
    xref.deleted = 'N' AND
    xref.ac = ref_map.accession AND
    ref_map.reference_id = ref.id AND
    {terms}
"""

SINGLE_SQL = BASE_SQL.format(
    terms="xref.upi = %(upi)s AND xref.taxid = %(taxid)s"
)

RANGE_SQL = BASE_SQL.format(terms="rna.id BETWEEN %(min_id)s AND %(max_id)s")


def fetch_dicts(cursor):
    """
    Creates a generator of all values as dicts for the results of the given
    query.
    """

    names = [col.name for col in cursor.description]
    while True:
        results = cursor.fetchmany()
        if not results:
            break
        for result in results:
            yield dict(zip(names, result))


def export_range(cursor, min_id, max_id):
    """
    Generates a series of XML strings representing all entries in the given
    range of ids.
    """

    with cursor() as cur:
        cur.execute(RANGE_SQL, {'min_id': min_id, 'max_id': max_id})
        results = fetch_dicts(cursor)
        grouped = it.groupby(results, op.itemgetter('upi', 'taxid'))
        for _, results in grouped:
            yield XmlEntry.build(list(results))


def export_upi(cursor, upi, taxid):
    """
    Will create a XmlEntry object for the given upi, taxid.
    """

    with cursor() as cur:
        cur.execute(SINGLE_SQL, {'upi': upi, 'taxid': taxid})
        data = list(fetch_dicts(cur))
        return XmlEntry.build(data)


def as_document(results):
    """
    This will create the required root XML element and place all the given
    XmlEntry objects as ElementTree.Element's in it. This then produces the
    string representation of that document which can be saved.
    """

    root = ET.Element('database')
    ET.SubElement(root, 'name', text='RNAcentral')
    ET.SubElement(
        root,
        'description',
        text='a database for non-protein coding RNA sequences'
    )
    ET.SubElement(root, 'release', text='1.0')

    timestamp = date.today.strftime('%d/%m/%Y')
    ET.SubElement(root, 'release_date', timestamp)

    count = 0
    count_element = ET.SubElement(root, 'entry_count')
    entries = ET.SubElement(root, 'entries')
    for result in results:
        entries.append(result.as_xml())
        count += 1

    if not count:
        raise ValueError("No entries found")

    count_element.text = count
    return ET.tostring(root)
