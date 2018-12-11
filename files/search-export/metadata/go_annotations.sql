COPY (
SELECT
json_build_object(
  'rna_id', anno.rna_id,
  'go_annotations', array_agg(
    json_build_object(
      'go_term_id', anno.ontology_term_id,
      'qualifier', anno.qualifier,
      'go_name', ont.name,
      'assigned_by', anno.assigned_by
    )
  )
)
FROM go_term_annotations anno
JOIN ontology_terms ont ON ont.ontology_term_id = anno.ontology_term_id
GROUP BY anno.rna_id
) TO STDOUT
