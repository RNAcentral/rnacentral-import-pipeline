BEGIN TRANSACTION;

INSERT INTO rnc_locus (
  taxid,
  assembly_id,
  chromosome,
  strand,
  locus_name,
  locus_start,
  locus_stop
) ( 
SELECT
  taxid,
  assembly_id,
  chromosome,
  strand,
  locus_start,
  locus_stop
FROM load_locus
) ON CONFLICT (locus_name) DO NOTHING;

-- We remove all existing assignments because it is possible that a sequence
-- will get moved from one locus to another. This 
DELETE 
FROM rnc_locus_members members
JOIN load_locus load
ON
  load.urs_taxid = members.urs_taxid
;

INSERT INTO rnc_locus_members (
  urs_taxid,
  region_id,
  locus_id,
  is_representative
) ( 
SELECT
  load.urs_taxid,
  load.region_id,
  locus.id,
  load.is_representative
FROM load_genes load
JOIN rnc_locus locus
ON
  locus.locus_name = load.locus_name
);

UPDATE rnc_rna_precomputed
SET
  is_locus_representative = bool_or(members.is_representative)
FROM rnc_locus_members

COMMIT;
