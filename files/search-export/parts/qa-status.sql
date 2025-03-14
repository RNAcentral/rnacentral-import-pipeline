COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'has_issue', qa.has_issue,
      'possible_contamination', qa.possible_contamination,
      'incomplete_sequence', qa.incomplete_sequence,
      'missing_rfam_match', qa.missing_rfam_match,
      'possible_orf', qa.possible_orf
    )
  FROM search_export_urs todo
  JOIN qa_status qa ON qa.rna_id = todo.urs_taxid
  ORDER BY todo.id
) TO STDOUT
