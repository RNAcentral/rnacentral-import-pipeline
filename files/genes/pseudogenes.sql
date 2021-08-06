COPY (
SELECT
  json_build_object(
      'assembly_id', :'assembly_id',
      'region_id', max(regions.region_name),
      'rna_id', max(regions.gene),
      'rna_type',  'pseudogene',
      'databases', ARRAY['ensembl'],
      'providing_databases', ARRAY['ensembl'],
      'chromosome', max(regions.chromosome),
      'strand', max(regions.strand),
      'identity', null,
      'was_mapped', false,
      'exons', array_agg(distinct exons.*)
  )
FROM ensembl_pseudogene_regions regions
JOIN ensembl_pseudogene_exons exons
ON
  exons.region_id = regions.id
WHERE
  regions.assembly_id = :'assembly_id'
GROUP BY regions.id
ORDER BY max(regions.chromosome), regions.region_start, regions.id
) TO STDOUT
