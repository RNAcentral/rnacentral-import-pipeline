COPY (
select
  json_build_object(
    'id', todo.id,
    'urs_id', todo.precompute_urs_id,
    'urs_taxid', todo.urs_taxid,
    'assembly_id', region.assembly_id,
    'chromosome', region.chromosome,
    'strand', region.strand,
    'start', region.region_start,
    'stop', region.region_stop
  )
FROM precompute_urs_taxid todo
JOIN rnc_sequence_regions region 
ON 
  region.urs_taxid = todo.urs_taxid
order by todo.precompute_urs_id, todo.id
) TO STDOUT
