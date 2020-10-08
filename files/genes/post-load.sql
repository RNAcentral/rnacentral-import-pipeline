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
SELECT DSTINCT
  taxid,
  assembly_id,
  chromosome,
  strand,
  locus_name,
  locus_start,
  locus_stop
FROM load_locus
) ON CONFLICT (assembly_id, locus_name) DO NOTHING;

-- We remove all existing assignments because it is possible that a sequence
-- will get moved from one locus to another. This 
DELETE 
FROM rnc_locus_members members
USING load_locus load
WHERE
  load.urs_taxid = members.urs_taxid
;

INSERT INTO rnc_locus_members (
  urs_taxid,
  region_id,
  locus_id,
  is_representative
) ( 
SELECT DISTINCT
  load.urs_taxid,
  load.region_id,
  locus.id,
  load.is_representative
FROM load_locus load
JOIN rnc_locus locus
ON
  locus.locus_name = load.locus_name
);

UPDATE rnc_rna_precomputed pre
SET
  is_locus_representative = t.is_representative
FROM (
  SELECT 
    urs_taxid urs_taxid,
    bool_or(is_representative) is_representative
  FROM rnc_locus_members members
  GROUP BY members.urs_taxid
) t
WHERE
  t.urs_taxid = pre.id
;

DROP TABLE load_locus;

COMMIT;
