COPY (
SELECT
json_build_object(
  'id', anno.rna_id,
  'go_term_id', anno.ontology_term_id,
  'qualifier', anno.qualifier,
  'go_name', ont.name,
  'assigned_by', anno.assigned_by
)
FROM go_term_annotations anno
JOIN ontology_terms ont ON ont.ontology_term_id = anno.ontology_term_id
) TO STDOUT
