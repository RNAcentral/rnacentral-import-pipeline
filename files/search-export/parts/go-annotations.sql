COPY (
SELECT
json_build_object(
  'id', todo.id,
  'urs_taxid', todo.urs_taxid,
  'go_term_id', anno.ontology_term_id,
  'qualifier', anno.qualifier,
  'go_name', ont.name,
  'assigned_by', anno.assigned_by
)
FROM search_export_urs todo
JOIN go_term_annotations anno
ON
  anno.rna_id = todo.urs_taxid
JOIN ontology_terms ont 
ON ont.ontology_term_id = anno.ontology_term_id
ORDER BY todo.id
) TO STDOUT
