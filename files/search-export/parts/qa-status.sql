COPY (
  SELECT
    json_build_object(
      'urs_taxid', qa.urs_taxid,
      'has_issue', qa.has_issue,
      'possible_contamination', qa.possible_contamination,
      'incomplete_sequence', qa.incomplete_sequence,
      'missing_rfam_match', qa.missing_rfam_match
    )
  FROM search_export_urs todo
  JOIN qa_status qa ON qa.urs_taxid = todo.urs_taxid
  ORDER BY todo.id
) TO STDOUT
