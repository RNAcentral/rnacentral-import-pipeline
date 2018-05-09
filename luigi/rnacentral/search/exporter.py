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

import json
import subprocess as sp
from datetime import date

from lxml import etree
from lxml.builder import E

from rnacentral.psql import PsqlWrapper

from .data import builder

XML_SCHEMA = 'http://www.ebi.ac.uk/ebisearch/XML4dbDumps.xsd'


BASE_SQL = """
SELECT
    json_build_object(
        'upi', rna.upi,
        'taxid', xref.taxid,
        'first_seen', array_agg(release1.timestamp),
        'last_seen', array_agg(release2.timestamp),
        'cross_references', array_agg(
            json_build_object(
                'name', acc."database",
                'external_id', acc.external_id,
                'optional_id', acc.optional_id,
                'accession', acc.accession,
                'non_coding_id', acc.non_coding_id,
                'parent_accession', acc.parent_ac || '.' || acc.seq_version
                )
            ),
        'description', array_agg(pre.description),
        'deleted', array_agg(xref.deleted),
        'length', array_agg(rna.len),
        'species', array_agg(acc.species),
        'organelles', array_agg(acc.organelle),
        'expert_dbs', array_agg(db.display_name),
        'rna_type', array_agg(pre.rna_type),
        'product', array_agg(acc.product),
        'md5', array_agg(rna.md5),
        'authors', array_agg(refs.authors),
        'journals', array_agg(refs.location),
        'pub_titles', array_agg(refs.title),
        'pub_ids', array_agg(refs.id),
        'pubmed_ids', array_agg(pubmed.ref_pubmed_id::varchar) || array_agg(refs.pmid),
        'dois', array_agg(pubmed.doi) || array_agg(refs.doi),
        'has_coordinates', array_agg(pre.has_coordinates),
        'rfam_family_names', array_agg(models.short_name),
        'rfam_ids', array_agg(hits.rfam_model_id),
        'rfam_clans', array_agg(models.rfam_clan_id),
        'rfam_status',
            case
                when cardinality((array_agg(pre.rfam_problems))) = 0 then '{{}}'
                else (array_agg(pre.rfam_problems))[1]::json
            end,
        'tax_strings', array_agg(acc.classification),
        'functions', array_agg(acc.function),
        'genes', array_agg(acc.gene),
        'gene_synonyms', array_agg(acc.gene_synonym),
        'common_name', array_agg(acc.common_name),
        'notes', array_agg(acc.note),
        'locus_tags', array_agg(acc.locus_tag),
        'standard_names', array_agg(acc.standard_name),
        'products', array_agg(acc.product),
        'go_annotations', array_agg(
            json_build_object(
                'go_term_id', anno.ontology_term_id,
                'qualifier', anno.qualifier,
                'go_name', ont.name,
                'assigned_by', anno.assigned_by
            )
        )
    )
FROM xref xref
JOIN rnc_accessions acc ON xref.ac = acc.accession
JOIN rnc_database db ON xref.dbid = db.id
JOIN rnc_release release1 ON xref.created = release1.id
JOIN rnc_release release2 ON xref.last = release2.id
JOIN rna rna ON xref.upi = rna.upi
JOIN rnc_rna_precomputed pre
ON
    xref.upi = pre.upi
    AND xref.taxid = pre.taxid
LEFT JOIN rnc_reference_map ref_map ON ref_map.accession = acc.accession
LEFT JOIN rnc_references refs ON refs.id = ref_map.reference_id
LEFT JOIN rfam_model_hits hits ON xref.upi = hits.upi
LEFT JOIN rfam_models models
ON
    hits.rfam_model_id = models.rfam_model_id
LEFT JOIN go_term_annotations anno ON anno.rna_id = pre.id
LEFT JOIN go_term_publication_map go_map
ON
    go_map.go_term_annotation_id = anno.go_term_annotation_id
LEFT JOIN ref_pubmed pubmed ON pubmed.ref_pubmed_id = go_map.ref_pubmed_id
LEFT JOIN ontology_terms ont
ON
    ont.ontology_term_id = anno.ontology_term_id
WHERE
  xref.deleted = 'N'
  AND %s
GROUP BY rna.upi, xref.taxid
"""

SINGLE_SQL = BASE_SQL % "xref.upi = '{upi}' AND xref.taxid = {taxid}"

RANGE_SQL = BASE_SQL % "rna.id BETWEEN {min_id} AND {max_id}"


def export(db, query, **kwargs):
    psql = PsqlWrapper(db)
    for result in psql.copy_to_iterable(query, **kwargs):
        try:
            data = json.loads(result['json_build_object'])
            yield builder(data)
        except:
            raise


def range(db, min_id, max_id):
    """
    Generates a series of XML strings representing all entries in the given
    range of ids.
    """
    return export(db, RANGE_SQL, min_id=min_id, max_id=max_id)


def upi(db, upi, taxid):
    """
    Will create a XmlEntry object for the given upi, taxid.
    """

    results = export(db, SINGLE_SQL, upi=upi, taxid=taxid)
    try:
        return next(results)
    except StopIteration:
        raise ValueError("Found no entries for %s_%i" % (upi, taxid))


def write(handle, results):
    """
    This will create the required root XML element and place all the given
    XmlEntry objects as ElementTree.Element's in it. This then produces the
    string representation of that document which can be saved.
    """

    handle.write('<database>')
    handle.write(etree.tostring(E.name('RNAcentral')))
    handle.write(etree.tostring(E.description('a database for non-protein coding RNA sequences')))
    handle.write(etree.tostring(E.release('1.0')))
    handle.write(etree.tostring(E.release_date(date.today().strftime('%d/%m/%Y'))))

    count = 0
    handle.write('<entries>')
    for result in results:
        count += 1
        handle.write(etree.tostring(result))
    handle.write('</entries>')

    if not count:
        raise ValueError("No entries found")

    handle.write(etree.tostring(E.entry_count(str(count))))
    handle.write('</database>')


def validate(filename):
    """
    Run xmllint validation on the given filename.
    """

    cmd = ('xmllint', filename, '--schema', XML_SCHEMA, '--stream')
    sp.check_call(cmd)
