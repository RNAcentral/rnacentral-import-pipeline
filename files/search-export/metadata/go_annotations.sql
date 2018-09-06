COPY (
SELECT
json_build_object(
  'rna_id', anno.rna_id,
  'go_annotations', array_agg(
    json_build_object(
      'go_term_id', anno.ontology_term_id,
      'qualifier', anno.qualifier,
      'go_name', ont.name,
      'assigned_by', anno.assigned_by,
      'pubmed_ids', pubmed.ref_pubmed_id,
      'dois', pubmed.doi
    )
  )
)
FROM go_term_annotations anno
JOIN ontology_terms ont ON ont.ontology_term_id = anno.ontology_term_id
LEFT JOIN go_term_publication_map go_map
ON
  go_map.go_term_annotation_id = anno.go_term_annotation_id
LEFT JOIN ref_pubmed pubmed ON pubmed.ref_pubmed_id = go_map.ref_pubmed_id
GROUP BY anno.rna_id
) TO STDOUT
