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
  locus_public_name,
  chromosome,
  strand,
  locus_start,
  locus_stop,
  member_count,
) ( 
SELECT DSTINCT
  assembly_id,
  locus_name,
  locus_public_name,
  chromosome,
  strand,
  locus_start,
  locus_stop,
  member_count,
FROM load_locus
) ON CONFLICT (assembly_id, locus_name) DO NOTHING;

INSERT INTO rnc_locus_members (
  urs_taxid,
  region_id,
  locus_id,
  member_type
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

DROP TABLE load_locus;

INSERT INTO rnc_gene_status (
  urs_taxid,
  region_id,
  assembly_id,
  cluster_status
) (
SELECT
  load.urs_taxid,
  load.region_id,
  load.assembly_id,
  load.cluster_status
FROM load_gene_status
);

COMMIT;
