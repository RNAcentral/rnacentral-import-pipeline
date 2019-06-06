COPY (
  SELECT 
    rfam_model_id 
  FROM rfam_models
  WHERE
    rfam_rna_type != 'Gene; rRNA'
) TO STDOUT;
