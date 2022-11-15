# -*- coding: utf-8 -*-

"""
Copyright [2009-current] EMBL-European Bioinformatics Institute
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
from rnacentral_pipeline import utils
from rnacentral_pipeline.rnacentral import lookup

import psycopg2
import psycopg2.extras

QUERY = """
select
	pre.id as urs_taxid,
	pre.rna_type,
	ra.rna_type so_term,
	COALESCE(rna.seq_short, rna.seq_long) as sequence,
	pre.description,
	ra.optional_id,
	external_id
from xref join
rnc_accessions ra on
ra.accession = xref.ac
join rna
on rna.upi = xref.upi
join rnc_rna_precomputed pre
on rna.upi = pre.upi
where ra.optional_id in %(aliases)s
or external_id in %(aliases)s
"""

EXON_QUERY = """
select exon_start, exon_stop, strand, chromosome from rnc_sequence_regions sr
join rnc_sequence_exons ex on ex.region_id = sr.id
where urs_taxid = %(urs_taxid)s
"""

def mapping(db_url, data):
	"""
	data is a dict of ext_id : aliases

	Want to get at least one hit from the query
	"""
	mapping = as_mapping(db_url, data)
	for value in mapping.values():
		value["sequence"] = value["sequence"].replace("U", "T")
	return mapping

def as_mapping(db_url, data):
	mapping = {}
	for ext_id in data:
		aliases = [a.strip() for a in data[ext_id] if a.strip() != '']
		for result in run_query(db_url, {'aliases':tuple(aliases),}, QUERY):
			exons = list(run_query(db_url, {'urs_taxid':result['urs_taxid'],}, EXON_QUERY))
			result['exons'] = exons
			mapping[ext_id] = result
			if '_' in result['urs_taxid']:
				break
	return mapping


def run_query(db_url, substitutions, query):
	data = {}
	conn = psycopg2.connect(db_url)

	cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
	cur.execute(query, substitutions)
	for result in cur:
		yield result
	cur.close()
	conn.close()
