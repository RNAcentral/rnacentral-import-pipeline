DROP TABLE IF EXISTS :tablename;

CREATE TABLE :tablename AS
SELECT 
	t.upi,
	models.rfam_model_id as model
FROM (
SELECT
	hits.upi upi,
	min(hits.rfam_hit_id) rfam_hit_id,
	min(hits.rfam_model_id) rfam_model_id
FROM rfam_model_hits hits
GROUP BY hits.upi
HAVING count(hits.rfam_hit_id) = 1
) t
JOIN rfam_models models ON t.rfam_model_id = models.rfam_model_id
WHERE
  rfam_rna_type != 'Gene; rRNA'
;

CREATE UNIQUE INDEX un_traveler_sequences_to_analyze__upi ON :tablename(upi);

INSERT INTO :tablename (
  upi,
  model
) (
SELECT DISTINCT
  pre.upi,
  'rRNA' as model
FROM rnc_rna_precomputed pre
JOIN qa_status qa ON qa.rna_id = pre.id
WHERE
  pre.is_active = true
  AND pre.rna_type = 'rRNA'
  AND qa.incomplete_sequence = false
) ON CONFLICT (upi) DO UPDATE
SET
  model = excluded.model
;

CREATE INDEX ix_traveler_sequences_to_analyze__model ON :tablename(model);

DELETE FROM :tablename urs
USING rnc_secondary_structure_layout layout
WHERE
  layout.urs = urs.upi
;
