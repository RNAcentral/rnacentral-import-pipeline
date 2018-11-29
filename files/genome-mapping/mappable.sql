COPY (
SELECT
  ensembl_url,
  assembly_id,
  taxid,
  division
FROM ensembl_assembly
WHERE
  division != 'EnsemblProtists'
  AND division != 'EnsemblFungi'
) TO STDOUT CSV;
