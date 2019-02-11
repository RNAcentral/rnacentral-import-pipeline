ALTER TABLE load_compara
  ADD COLUMN urs_taxid text,
  ADD COLUMN homology_id int
;

-- Determine all the urs_taxids to store
UPDATE load_compara
SET
  urs_taxid = xref.upi || '_' || xref.taxid
FROM xref
WHERE
  xref.accession = load_compara.ensembl_transcript
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

-- Remove all old compara data so we don't have stale data.
TRUNCATE TABLE ensembl_compara;

INSERT INTO ensembl_compara (
  urs_taxid,
  ensembl_transcript_id,
  homology_group
) (
SELECT DISTINCT
  load.urs_taxid,
  load.ensembl_transcript,
  load.homology_group
FROM load_compara load
)
;

DROP TABLE load_compara;
