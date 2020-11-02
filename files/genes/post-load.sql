BEGIN TRANSACTION;

-- DELETE all assignments from the assemblies to be loaded. This is needed
-- because
DELETE
FROM rnc_locus locus
USING load_locus load
WHERE
  load.assembly_id = locus.assembly_id
;

DELETE
FROM rnc_gene_status status
USING load_gene_status load
WHERE
  load.urs_taxid = status.urs_taxid
;

INSERT INTO rnc_locus (
  assembly_id,
  locus_name,
  public_locus_name,
  chromosome,
  strand,
  locus_start,
  locus_stop,
  member_count
) ( 
SELECT DISTINCT
  assembly_id,
  locus_name,
  locus_public_name,
  chromosome,
  strand,
  locus_start,
  locus_stop,
  member_count,
FROM load_locu
) ON CONFLICT (assembly_id, locus_name) DO NOTHING;

INSERT INTO rnc_locus_members (
  urs_taxid,
  region_id,
  locus_id,
  membership_status
) ( 
SELECT DISTINCT
  load.urs_taxid,
  load.region_id,
  locus.id,
  load.member_type
FROM load_locus load
JOIN rnc_locus locus
ON
  locus.locus_name = load.locus_name
);

INSERT INTO rnc_gene_status (
  urs_taxid,
  region_id,
  assembly_id,
  status
) (
SELECT
  load.urs_taxid,
  load.region_id,
  load.assembly_id,
  load.status
FROM load_gene_status load
);

INSERT INTO rnc_gene_status (
  assembly_id, 
  urs_taxid, 
  region_id, 
  status
) (
SELECT 
  locus.assembly_id, 
  member.urs_taxid, 
  member.region_id, 
  'clustered' 
FROM rnc_locus locus 
JOIN rnc_locus_members member 
ON member.locus_id = locus.id
);

DROP TABLE load_locus;
DROP TABLE load_gene_status;

COMMIT;
