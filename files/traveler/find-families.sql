COPY (
  SELECT distinct
    rfam_model_id 
  FROM urs_with_one_rfam
  WHERE
    rfam_rna_type != 'Gene; rRNA'
) TO STDOUT
