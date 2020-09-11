COPY (
select
  json_build_object(
    'id', todo.id,
    'chromosomes', ARRAY(
      select
        distinct chromosome
      from rnc_sequence_regions
      where urs_taxid = todo.urs_taxid
    ),
    'has_coordinates', exists(
        select 1
        from rnc_sequence_regions
        where urs_taxid = todo.urs_taxid
    ),
  )
FROM :tablename todo
) TO STDOUT
