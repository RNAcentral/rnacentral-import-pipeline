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

import operator as op
import itertools as it
import xml.etree.ElementTree as ET

from .data import XmlEntry


SQL = """
SELECT
    xref.taxid,
    xref.deleted,
    acc.species,
    acc.organelle,
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
    acc.locus_tag,
    acc.standard_name,
    ref.id pub_id,
    ref.location,
    ref.title,
    ref.pmid,
    ref.doi,
    ref.author
FROM
    xref,
    rnc_accessions acc,
    rnc_database db,
    rnc_release release1,
    rnc_release release2,
    rna,
    rnc_rna_precomputed pre,
    rnc_references refs,
    rnc_reference_map ref_map
WHERE
    xref.ac = acc.accession AND
    xref.dbid = db.id AND
    xref.created = release1.id AND
    xref.last = release2.id AND
    xref.upi = rna.upi AND
    xref.upi = pre.upi AND
    xref.taxid = pre.taxid AND
    xref.upi = :upi AND
    xref.deleted = 'N' AND
    xref.ac = ref_map.accession AND
    ref_map.reference = ref.id AND
    xref.taxid = :taxid
"""


def export(cursor, min_id, max_id):
    """
    Generates a series of XML strings representing all entries in the given
    range of ids.
    """

    with cursor() as cur:
        results = cur.execute(SQL, min_id=min_id, max_id=max_id)
        grouped = it.groupby(results, op.attrgetter('upi', 'taxid'))
        for _, results in grouped:
            entry = XmlEntry.build(results)
            yield ET.tostring(entry.as_xml())
