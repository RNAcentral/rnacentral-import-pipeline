COPY (
  SELECT
    ensembl_url,
    assembly_id,
    taxid,
    gca_accession
  FROM ensembl_assembly
  WHERE
    selected_genome = true
) TO STDOUT CSV;
