COPY (
  SELECT
    ensembl_url,
    assembly_id,
    taxid,
    division
  FROM ensembl_assembly
  WHERE
    selected_genome = true
) TO STDOUT CSV;
