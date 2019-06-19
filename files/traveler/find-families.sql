DROP TABLE IF EXISTS urs_with_one_rfam;

CREATE TABLE urs_with_one_rfam AS
SELECT 
	t.upi,
	t.rfam_hit_id,
	models.*
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
;

CREATE INDEX ix_urs_with_one_rfam__upi ON urs_with_one_rfam(upi);
CREATE INDEX ix_urs_with_one_rfam__rfam_hit_id ON urs_with_one_rfam(rfam_hit_id);

COPY (
  SELECT 
    rfam_model_id 
  FROM rfam_models
  WHERE
    rfam_rna_type != 'Gene; rRNA'
) TO STDOUT;
