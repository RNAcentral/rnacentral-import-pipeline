COPY (
SELECT
  species.ensembl_url,
  species.assembly_id,
  species.division
FROM ensembl_assembly
) TO STDOUT CSV;

