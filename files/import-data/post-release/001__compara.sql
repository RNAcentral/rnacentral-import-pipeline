ALTER TABLE load_compara
  ADD COLUMN urs_taxid text
;

UPDATE load_compara
SET
  urs_taxid = xref.upi || '_' || xref.taxid
FROM xref
WHERE
  xref.accession = load_compara.ensembl_transcript
;

DELETE FROM load_compara 
WHERE urs_taxid IS NULL
;

DELETE FROM ensembl_compara compara
USING load_compara load
WHERE
  load.urs_taxid = compara.urs_taxid
;

INSERT INTO ensembl_compara (
  homology_group,
  urs_taxid
) (
SELECT DISTINCT
  load.homology_group,
  load.urs_taxid
FROM load_compara load
)
;

DROP TABLE load_compara;
