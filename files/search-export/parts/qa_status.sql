COPY (
  SELECT
    json_build_object(
      'id', qa.urs_taxid,
      'has_issue', qa.has_issue,
      'possible_contamination', qa.possible_contamination,
      'incomplete_sequence', qa.incomplete_sequence,
      'missing_rfam_match', qa.missing_rfam_match
    )
  FROM qa_status qa
) TO STDOUT
