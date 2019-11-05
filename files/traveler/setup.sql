DROP TABLE IF EXISTS :tablename;

CREATE TABLE :tablename AS
SELECT 
	t.upi,
  case
  when clans.rfam_clan_id == 'CL00111' THEN 'crw'
  when clans.rfam_clan_id == 'CL00112' and models.rfam_model_id != 'CL00112' THEN 'ribovision'
  else models.rfam_model_id
  end as model
FROM (
SELECT
	hits.upi upi,
	min(hits.rfam_hit_id) rfam_hit_id,
	min(hits.rfam_model_id) rfam_model_id
FROM rfam_model_hits hits
GROUP BY hits.upi
HAVING count(hits.rfam_hit_id) = 1
) t
JOIN rfam_models models 
ON 
  t.rfam_model_id = models.rfam_model_id
LEFT JOIN rfam_clans clans
ON
  models.rfam_clan_id = clans.rfam_clan_id
;

CREATE UNIQUE INDEX un_traveler_sequences_to_analyze__upi ON :tablename(upi);
CREATE INDEX ix_traveler_sequences_to_analyze__model ON :tablename(model);

-- Delete already computed
DELETE FROM :tablename urs
USING rnc_secondary_structure_layout layout
WHERE
  layout.urs = urs.upi
;

DELETE FROM :tablename to_draw
JOIN qa_status qa 
ON 
  qa.upi = to_draw.upi
WHERE
  to_draw.model in ('crw', 'ribovision')
  AND qa.incomplete_sequence = true
;

-- Delete all lncRNAs
DELETE FROM :tablename urs
USING rnc_rna_precomputed pre
WHERE
  layout.urs = pre.upi
  and pre.rna_type = 'lncRNA'
;
