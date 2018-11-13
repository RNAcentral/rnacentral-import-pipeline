COPY (
SELECT
  ensembl_url,
  assembly_id,
  taxid,
  division
from ensembl_assembly
where
  common_name = 'human'
) TO STDOUT CSV;
