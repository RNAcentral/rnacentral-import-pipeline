COPY (
  select
  json_build_object(
    'taxid', max(pre.taxid),
    'urs_taxid', max(pre.id),
    'assembly_id', :'assembly_id',
    'region_id', max(regions.id),
    'region_name', max(regions.region_name),
    'region_start', max(regions.region_start),
    'region_stop', max(regions.region_stop),
    'insdc_rna_type',  max(pre.rna_type),
    'so_rna_type', max(ont.name),
    'databases', regexp_split_to_array(max(pre."databases"), ','),
    'providing_databases', max(regions.providing_databases),
    'chromosome', max(regions.chromosome),
    'strand', max(regions.strand),
    'identity', max(regions.identity),
    'was_mapped', bool_or(regions.was_mapped),
    'exons', array_agg(distinct exons.*),
    'qa', array_agg(json_build_object(
      'has_issue', qa.has_issue,
      'incomplete_sequence', qa.incomplete_sequence,
      'possible_contamination', qa.possible_contamination,
      'missing_rfam_match', qa.missing_rfam_match
    ))
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
JOIN ontology_terms ont
ON
  ont.ontology_term_id = pre.so_rna_type
WHERE
  pre.is_active = true
  AND regions.assembly_id = :'assembly_id'
GROUP BY regions.id
ORDER BY max(regions.chromosome), max(regions.strand), max(pre.rna_type)
) TO STDOUT
