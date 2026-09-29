\timing

BEGIN;

ALTER TABLE load_compara
  ADD COLUMN urs_taxid text,
  ADD COLUMN homology_id int,
  ADD COLUMN dbid int
;

-- Ensembl ncRNA (25) and Ensembl mRNA (61) transcripts, retired ones included
-- so their old rows can be replaced. The mRNA database keeps each transcript's
-- UTRs under the same external_id, so only its mRNA entry stands for it.
CREATE TEMP TABLE compara_transcripts AS
select
  xref.urs || '_' || xref.taxid as urs_taxid,
  acc.external_id as transcript,
  xref.dbid,
  xref.deleted
from xref
join rnc_accessions acc on acc.accession = xref.ac
where
  xref.dbid IN (25, 61)
  and (xref.dbid = 25 or acc.feature_name = 'mRNA')
;

-- Determine all the urs_taxids to store
UPDATE load_compara
SET
  urs_taxid = t.urs_taxid,
  dbid = t.dbid
FROM compara_transcripts t
WHERE
  t.deleted = 'N'
  and t.transcript = load_compara.ensembl_transcript
;

-- populate the load table with the required homology ids.
UPDATE load_compara load
SET
  homology_id = t.homology_id
FROM (
  select
    homology_group,
    nextval('ensembl_compara_homology_id') as homology_id
  from load_compara
  group by homology_group
) as t
where
  t.homology_group = load.homology_group
;

DELETE FROM load_compara
WHERE urs_taxid IS NULL
;

-- Replace only the databases in this load, so loading mRNA homologies keeps
-- the ncRNA ones and the other way round.
DELETE FROM ensembl_compara compara
USING compara_transcripts t
WHERE
  t.transcript = compara.ensembl_transcript_id
  and t.dbid IN (SELECT DISTINCT dbid FROM load_compara)
;

-- Drop indexes before bulk insert; both are recreated below
DROP INDEX IF EXISTS rnacen.fk_ensembl_compara__urs_taxid;
DROP INDEX IF EXISTS rnacen.ix_ensembl_compara__homology_id;

INSERT INTO ensembl_compara (
  urs_taxid,
  ensembl_transcript_id,
  homology_id
) (
SELECT DISTINCT
  load.urs_taxid,
  load.ensembl_transcript,
  load.homology_id
FROM load_compara load
)
;

DROP TABLE load_compara;

-- Recreate indexes
CREATE INDEX fk_ensembl_compara__urs_taxid ON rnacen.ensembl_compara USING btree (urs_taxid);
CREATE INDEX ix_ensembl_compara__homology_id ON rnacen.ensembl_compara USING btree (homology_id);

COMMIT;
