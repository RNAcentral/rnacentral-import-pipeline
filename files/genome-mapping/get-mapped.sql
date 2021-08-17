COPY (
  SELECT DISTINCT
    urs_taxid
  FROM pipeline_tracking_genome_mapping
  WHERE
    assembly_id = :'assembly_id'
) TO STDOUT CSV;
