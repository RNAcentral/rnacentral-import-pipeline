COPY (
  SELECT DISTINCT
    urs_taxid
  FROM rnc_pipeline_tracking
  WHERE
    assembly_id = :'assembly_id'
) TO STDOUT CSV;
