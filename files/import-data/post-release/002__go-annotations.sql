INSERT INTO go_term_annotations (
    rna_id,
    qualifier,
    ontology_term_id,
    evidence_code,
    assigned_by,
    extensions
) (
SELECT
    rna_id,
    qualifier,
    ontology_term_id,
    evidence_code,
    assigned_by,
    extensions
FROM load_go_term_annotations
)
ON CONFLICT (rna_id, qualifier, ontology_term_id, evidence_code, assigned_by)
DO UPDATE
SET
    extensions = excluded.extensions
;

DROP TABLE load_go_term_annotations;
