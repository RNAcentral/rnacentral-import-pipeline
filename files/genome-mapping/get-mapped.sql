COPY (
  SELECT DISTINCT
    urs_taxid
  FROM pipeline_tracking_genome_mapping
  WHERE
    assembly_id = :'assembly_id'
  UNION
    SELECT DISTINCT
      urs_taxid
    FROM rnc_sequence_regions
    WHERE
      assembly_id = :'assembly_id'
      AND was_mapped = false
) TO STDOUT CSV;
