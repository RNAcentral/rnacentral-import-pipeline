COPY (
SELECT
  ensembl_url,
  assembly_id,
  taxid,
  division
from ensembl_assembly
) TO STDOUT CSV;
