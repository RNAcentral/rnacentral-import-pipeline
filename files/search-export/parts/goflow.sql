COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'should_show_goflow', true
    )
    FROM search_export_urs todo
    JOIN go_flow_llm_curation_results gfllm
    ON
      todo.urs_taxid = gfllm.urs_taxid
    ORDER by todo.id
) TO STDOUT
