COPY (
  select
  json_build_object(
    'taxid', max(pre.taxid),
    'urs_taxid', max(pre.id),
    'assembly_id', :'assembly_id',
    'region_id', max(regions.id),
    'region_name', max(regions.region_name),
    'insdc_rna_type',  max(pre.rna_type),
    'so_rna_type', max(pre.so_rna_type),
    'databases', regexp_split_to_array(max(pre."databases"), ','),
    'providing_databases', max(regions.providing_databases),
    'chromosome', max(regions.chromosome),
    'strand', max(regions.strand),
    'identity', max(regions.identity),
    'was_mapped', bool_or(regions.was_mapped),
    'exons', array_agg(distinct exons.*),
    'qa', array_agg(distinct qa.*)
  ) 
FROM rnc_rna_precomputed pre
JOIN rnc_sequence_regions regions
ON
  regions.urs_taxid = pre.id
JOIN rnc_sequence_exons exons
ON
  exons.region_id = regions.id
JOIN qa_status qa
ON
  qa.rna_id = pre.id
WHERE
  pre.is_active = true
  AND regions.assembly_id = :'assembly_id'
GROUP BY regions.id
ORDER BY max(regions.chromosome), max(regions.strand), max(pre.rna_type)
) TO STDOUT
