INSERT INTO go_term_publication_map (
    go_term_annotation_id,
    ref_pubmed_id
) (
SELECT DISTINCT
    annotations.go_term_annotation_id,
    load_go_term_publication_map.pubmed_id
FROM load_go_term_publication_map
JOIN go_term_annotations annotations
ON
    annotations.rna_id = load_go_term_publication_map.rna_id
    AND annotations.qualifier = load_go_term_publication_map.qualifier
    AND annotations.assigned_by = load_go_term_publication_map.assigned_by
    AND annotations.ontology_term_id = load_go_term_publication_map.ontology_term_id
    AND annotations.evidence_code = load_go_term_publication_map.evidence_code
)
ON CONFLICT (go_term_annotation_id, ref_pubmed_id)
DO NOTHING
;

DROP TABLE load_go_term_publication_map;
