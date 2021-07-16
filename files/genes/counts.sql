COPY (
  select
  json_build_object(
    'urs_taxid', regions.urs_taxid,
    'assembly_id', :'assembly_id',
    'mapped_count', count(distinct regions.id) filter (where was_mapped = true),
    'given_count', count(distinct regions.id) filter (where was_mapped = false),
    'total_count', count(distinct regions.id)
  ) 
FROM rnc_rna_precomputed pre
JOIN rnc_sequence_regions regions
ON
  regions.urs_taxid = pre.id
WHERE
  pre.is_active = true
  AND regions.assembly_id = :'assembly_id'
  AND pre.taxid = :taxid
GROUP BY regions.urs_taxid, regions.assembly_id
) TO STDOUT
