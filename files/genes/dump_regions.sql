COPY(
  SELECT 
  sr.region_name as region_name, 
  sr.assembly_id as assembly_id,
  sr.id as region_id

  FROM rnc_sequence_regions sr 
  JOIN temp_gff_region_names rn ON rn.region_name = sr.region_name

  WHERE sr.urs_taxid LIKE '%_' || :'taxid'
) TO STDOUT WITH CSV HEADER;
