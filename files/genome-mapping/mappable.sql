COPY (
SELECT
  ensembl_url,
  assembly_id,
  url
from ensembl_assembly
) TO STDOUT CSV;
