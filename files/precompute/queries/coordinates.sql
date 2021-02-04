COPY (
select
  json_build_object(
    'id', todo.urs_taxid,
    'assembly_id', region.assembly_id,
    'chromosome', region.chromosome,
    'strand', region.strand,
    'start', region.region_start,
    'stop', region.region_stop
  )
FROM :tablename todo
JOIN rnc_sequence_regions region 
ON 
  region.urs_taxid = todo.urs_taxid
order by todo.id
) TO STDOUT
