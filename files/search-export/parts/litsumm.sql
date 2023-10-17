COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'should_show_litsumm', lss.should_show
    )
    FROM search_export_urs todo
    JOIN litsumm_summaries lss
    ON
      todo.urs_taxid = lss.primary_id
    ORDER by todo.id
) TO STDOUT
