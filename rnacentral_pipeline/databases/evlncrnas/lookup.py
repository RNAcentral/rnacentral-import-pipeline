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
import pandas as pd
import psycopg2
import psycopg2.extras

from rnacentral_pipeline import utils
from rnacentral_pipeline.rnacentral import lookup

QUERY = """
select DISTINCT ON (pre.id)
	pre.id as urs_taxid,
	pre.rna_type,
	ra.rna_type as so_term,
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
order by pre.id
"""

# create temp table ev_ids (
#     aliases TEXT
# );

# INSERT INTO ev_ids(aliases)
# VALUES( unnest(%(aliases)s ) );

STAGE_1_QUERY = """


create temp table ev_taxa (
    taxid INT
);

INSERT INTO ev_taxa(taxid)
VALUES( unnest(%(taxids)s ) );

with xref_cte(urs_taxid, ac) as (
	select
		xref.upi || '_' || xref.taxid as urs_taxid,
		xref.ac
	from xref join ev_taxa on xref.taxid = ev_taxa.taxid
	where xref.deleted = 'N'
),
ac_cte as (
	select
	accession,
	optional_id,
	external_id
	from rnc_accessions
	where external_id is not NULL
)

select
	xref_cte.urs_taxid,
	ac.optional_id,
	ac.external_id

from ac_cte ac
join
    xref_cte on ac.accession = xref_cte.ac
limit 10;

DROP TABLE ev_taxa;
"""

EXON_QUERY = """
select urs_taxid, exon_start, exon_stop, strand, chromosome from rnc_sequence_regions sr
join rnc_sequence_exons ex on ex.region_id = sr.id
join xref on xref.upi || '_' || xref.taxid = urs_taxid
where xref.deleted = 'N'
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
    # data = data.explode('Aliases').drop_duplicates(subset='Aliases').rename(columns={'Aliases':'external_id'})#.set_index('external_id')
    print(len(data))
    data = data.drop(data[data["Name"] == " "].index)
    print(data)
    taxa = (
        data["taxid"]
        .drop_duplicates()
        .sort_values()
        .dropna()
        .astype("int")
        .values.tolist()
    )
    aliases = data["Name"].drop_duplicates().sort_values().dropna().values.tolist()
    print(len(aliases), len(taxa))
    print(aliases[:10])
    db_mapping = run_query(db_url, {"taxids": taxa, "aliases": aliases}, STAGE_1_QUERY)
    print("Stage 1 query complete")
    # # urs_taxids = db_mapping['urs_taxid'].drop_duplicates().sort_values()
    # # print(len(urs_taxids))
    # # exons = run_query(db_url, {}, EXON_QUERY)
    # # print("Exon query complete")

    # # db_data_with_exons = db_mapping.join(exons, on="urs_taxid", how="inner")

    # print(db_mapping)

    exit()
    # print(data.columns, data)
    # print(db_mapping)
    # mapped_data = data.join(db_mapping, how='inner')

    # print(mapped_data)
    # return mapped_data


def run_query(db_url, substitutions, query):
    data = {}
    conn = psycopg2.connect(db_url)

    cur = conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor)
    cur.execute(query, substitutions)
    res = cur.fetchall()
    cur.close()
    conn.close()
    return pd.DataFrame(res)
