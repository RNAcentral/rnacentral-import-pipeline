BEGIN;

DELETE FROM rfam_go_terms rfam
USING load_rfam_go_terms load
WHERE
  load.rfam_model_id = rfam.rfam_model_id
;

INSERT INTO rfam_go_terms (
  ontology_term_id,
  rfam_model_id
) (
SELECT
  ontology_term_id,
  rfam_model_id
FROM load_rfam_go_terms
)
ON CONFLICT (ontology_term_id, rfam_model_id) DO UPDATE SET
  ontology_term_id = excluded.ontology_term_id,
  rfam_model_id = excluded.rfam_model_id
;

DROP TABLE load_rfam_go_terms;

COMMIT;
