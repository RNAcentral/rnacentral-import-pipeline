DROP TABLE IF EXISTS :tablename;
CREATE TABLE :tablename (
  urs_taxid text PRIMARY KEY,
  urs text NOT NULL,
  taxid int NOT NULL,
  seq text NOT NULL
);

INSERT INTO :tablename (
  urs,
  taxid,
  urs_taxid,
  seq
) (
SELECT
  rna.upi,
  pre.taxid,
  pre.id,
  COALESCE(rna.seq_short, rna.seq_long)
FROM rna
JOIN rnc_rna_precomputed pre 
ON 
  pre.upi = rna.upi
  and pre.taxid in (select taxid from ensembl_assembly)
WHERE 
  rna.len BETWEEN :min_length AND :max_length
);

DELETE FROM :tablename to_map
USING pipeline_tracking_genome_mapping gm
WHERE 
  gm.urs_taxid = to_map.urs_taxid
;
