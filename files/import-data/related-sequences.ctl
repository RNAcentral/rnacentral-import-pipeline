LOAD CSV
FROM stdin
HAVING FIELDS (
  source_accession,
  target_accession,
  relationship_type,
  methods
)
INTO {{PGDATABASE}}?load_rnc_related_sequences
TARGET COLUMNS (
  source_accession,
  target_accession,
  relationship_type,
  methods
)

WITH truncate,
  skip header = 0,
  fields escaped by double-quote,
  fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_rnc_related_sequences;
$$,
$$
create table load_rnc_related_sequences (
  source_accession varchar(100) NOT NULL,
  target_accession varchar(100) NOT NULL,
  relationship_type text NOT NULL,
  methods text[]
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_related_sequences (
  source_urs_taxid,
  source_accession,
  target_urs_taxid,
  target_accession,
  relationship_type,
  methods
) (
select
  source.id,
  load.source_accession,
  target.id,
  load.target_accession,
  load.relationship_type::related_sequence_relationship,
  load.methods
from load_rnc_related_sequences load
join xref source_xref on source_xref.ac = load.source_accession
join rnc_rna_precomputed source
ON
  source.upi = source_xref.upi
  and source.taxid = source_xref.taxid
left join xref target_xref on target_xref.ac = load.target_accession
left join rnc_rna_precomputed target
on
  target.upi = target_xref.upi
  and target.taxid = target_xref.taxid
)
ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
SET
  methods = EXCLUDED.methods || rnc_related_sequences.methods
;
$$,
$$
drop table load_rnc_related_sequences;
$$
;
