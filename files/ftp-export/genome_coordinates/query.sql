COPY (
SELECT
  json_build_object(
      'assembly_id', :'assembly_id',
      'region_id', max(regions.region_name),
      'rna_id', max(pre.id),
      'description', max(pre.short_description),
      'rna_type',  max(pre.rna_type),
      'databases', regexp_split_to_array(max(pre."databases"), ','),
      'providing_databases', array_agg(ac.database),
      'chromosome', max(regions.chromosome),
      'strand', max(regions.strand),
      'identity', max(regions.identity),
      'was_mapped', bool_or(regions.was_mapped),
      'exons', array_agg(distinct exons.*)
  )
FROM rnc_rna_precomputed pre
JOIN rnc_sequence_regions_active regions
ON
  regions.urs_taxid = pre.id
JOIN rnc_sequence_exons exons
ON
  exons.region_id = regions.id
LEFT JOIN rnc_accession_sequence_region sra on sra.region_id = regions.id
LEFT JOIN rnc_accessions ac on sra.accession = ac.accession
WHERE
  pre.is_active = true
  AND regions.assembly_id = :'assembly_id'
GROUP BY regions.id
ORDER BY max(regions.chromosome), max(regions.region_start), regions.id
) TO STDOUT
