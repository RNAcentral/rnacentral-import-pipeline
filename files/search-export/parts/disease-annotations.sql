COPY (
SELECT
json_build_object(
  'id', todo.id,
  'urs_taxid', todo.urs_taxid,
  'disease_term', anno.ontology_term_id,
  'disease_name', ont.name,
  'disease_text', anno.disease_text,
  'assigned_by', anno.assigned_by
)
FROM search_export_urs todo
JOIN rnc_disease_annotations anno
ON
  anno.urs_taxid = todo.urs_taxid
-- LEFT: an annotation whose text could not be mapped has no term, but still
-- counts as a disease annotation.
LEFT JOIN ontology_terms ont
ON ont.ontology_term_id = anno.ontology_term_id
ORDER BY todo.id
) TO STDOUT
