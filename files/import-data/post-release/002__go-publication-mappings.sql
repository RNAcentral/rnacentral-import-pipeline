INSERT INTO go_term_publication_map (
    go_term_annotation_id,
    reference_id
) (
SELECT DISTINCT
    annotations.go_term_annotation_id,
    refs.id
FROM load_go_term_publication_map load
JOIN go_term_annotations annotations
ON
    annotations.rna_id = load.rna_id
    AND annotations.qualifier = load.qualifier
    AND annotations.assigned_by = load.assigned_by
    AND annotations.ontology_term_id = load.ontology_term_id
    AND annotations.evidence_code = load.evidence_code
JOIN rnc_references refs
ON
  refs.pubmed_id = load.pubmed_id
)
ON CONFLICT (go_term_annotation_id, reference_id)
DO NOTHING
;

DROP TABLE load_go_term_publication_map;
