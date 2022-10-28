COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'publication_count', counts.publication_count
    )
    FROM search_export_urs todo
    JOIN search_export_publication_counts counts
    ON
      todo.urs_taxid = counts.urs
    ORDER by todo.id
) TO STDOUT
