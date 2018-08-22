COPY (
select
  assembly_id,
  ensembl_url,
  taxid
from ensembl_assembly
) TO STDOUT CSV
