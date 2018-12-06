-- Update the load table to contain the source URS
UPDATE load_rnc_related_sequences load
SET
  source_urs_taxid = xref.upi || '_' || xref.taxid
from xref
where
  xref.ac = load.source_accession
  and xref.deleted = 'N'
;

-- Remove things with no source URS, should do nothing
delete from load_rnc_related_sequences where source_urs_taxid is null;

-- Copy over all related proteins.
INSERT INTO rnc_related_sequences (
  source_urs_taxid,
  source_accession,
  target_accession,
  relationship_type,
  methods
) (
select
  load.source_urs_taxid,
  load.source_accession,
  load.target_accession,
  load.relationship_type::related_sequence_relationship,
  load.methods
from load_rnc_related_sequences load
WHERE
  load.relationship_type = 'target_protein'
)
ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
SET
  methods = EXCLUDED.methods || rnc_related_sequences.methods
;

-- Delete all target_proteins that should be copied over now
delete from load_rnc_related_sequences where relationship_type = 'target_protein';

-- For related RNA try to figure out what the accession/URS is.

-- If the accession is one we know just use it
INSERT INTO rnc_related_sequences (
  source_urs_taxid,
  source_accession,
  target_urs_taxid,
  target_accession,
  relationship_type,
  methods
) (
select distinct
  load.source_urs_taxid,
  load.source_accession,
  target.upi || '_' || target.taxid,
  load.target_accession,
  load.relationship_type::related_sequence_relationship,
  load.methods
from load_rnc_related_sequences load
join xref target on target.ac = load.target_accession
where
  load.relationship_type IN ('target_rna', 'isoform')
  and target.deleted = 'N'
)
ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
SET
  methods = EXCLUDED.methods || rnc_related_sequences.methods
;

-- Delete the loaded rnas with a known acccession
delete from load_rnc_related_sequences load
USING xref
WHERE
  xref.ac = load.target_accession
  and load.relationship_type = 'target_rna'
  and xref.deleted = 'N'
;

-- Build a table representing the related Ensembl genes
create temp table gene_upi_mapping as
select
  xref.upi || '_' || xref.taxid "urs_taxid",
  'ENSEMBL:' || split_part(acc.optional_id, '.', 1) "versionless_gene"
from xref
join rnc_accessions acc
ON
  acc.accession = xref.ac
where
  xref.deleted = 'N'
  and xref.dbid = 25
;

create index ix_gene_upi_mapping__versionless_gene on gene_upi_mapping(versionless_gene);

-- If the accession is a known ensembl gene copy that over
INSERT INTO rnc_related_sequences (
  source_urs_taxid,
  source_accession,
  target_urs_taxid,
  target_accession,
  relationship_type,
  methods
) (
select distinct on
  (load.source_accession, load.target_accession, load.relationship_type) load.source_urs_taxid,
  load.source_accession,
  gene.urs_taxid,
  load.target_accession,
  load.relationship_type::related_sequence_relationship,
  load.methods
from load_rnc_related_sequences load
join gene_upi_mapping gene
ON
  gene.versionless_gene = load.target_accession
where
  load.relationship_type = 'target_rna'
)
ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
SET
  methods = EXCLUDED.methods || rnc_related_sequences.methods
;

-- Cleanup the sequences with known gene
DELETE FROM load_rnc_related_sequences load
USING gene_upi_mapping gene
WHERE
  gene.versionless_gene = load.target_accession
  and load.relationship_type = 'target_rna'
;

-- Insert whatever remains with empty source_urs_taxid
INSERT INTO rnc_related_sequences (
  source_urs_taxid,
  source_accession,
  target_urs_taxid,
  target_accession,
  relationship_type,
  methods
) (
select distinct
  load.source_urs_taxid,
  load.source_accession,
  null,
  load.target_accession,
  load.relationship_type::related_sequence_relationship,
  load.methods
from load_rnc_related_sequences load
)
ON CONFLICT (source_accession, target_accession, relationship_type) DO UPDATE
SET
  methods = EXCLUDED.methods || rnc_related_sequences.methods
;

-- Ensure all methods are distinct
update rnc_related_sequences
set
  methods = ARRAY(select distinct unnest(methods))
;

drop table load_rnc_related_sequences;
