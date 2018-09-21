COPY (
SELECT
  json_build_object(
      'region_id', max(regions.region_name),
      'rna_id', max(pre.id),
      'rna_type',  max(pre.rna_type),
      'databases', regexp_split_to_array(max(pre."databases"), ','),
      'providing_databases', array_agg(distinct regions.providing_databases),
      'chromosome', max(regions.chromosome),
      'strand', max(regions.strand),
      'identity', max(regions.identity),
      'exons', array_agg(distinct exons.*)
  )
FROM rnc_rna_precomputed pre
JOIN rnc_sequence_regions regions
ON
  regions.urs_taxid = pre.id
JOIN rnc_sequence_exons exons
ON
  exons.region_id = regions.id
JOIN ensembl_coordinate_systems coords
ON
  coords.chromosome = regions.chromosome
WHERE
    pre.is_active = true
    AND regions.assembly_id = :assembly_id
    AND coords.is_reference = true
GROUP BY regions.id
ORDER BY max(coords.karyotype_rank), regions.id
) TO STDOUT
