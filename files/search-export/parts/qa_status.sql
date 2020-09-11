COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'has_issue', qa.has_issue,
      'possible_contamination', qa.possible_contamination,
      'incomplete_sequence', qa.incomplete_sequence,
      'missing_rfam_match', qa.missing_rfam_match
    )
  FROM :tablename todo
  JOIN qa_status qa
  ON
    qa.rna_id = todo.urs_taxid
) TO STDOUT
