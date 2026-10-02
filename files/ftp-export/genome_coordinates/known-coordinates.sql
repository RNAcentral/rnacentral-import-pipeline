COPY (
select distinct
  assembly.assembly_id,
  assembly.ensembl_url,
  assembly.taxid
from ensembl_assembly assembly
where assembly.selected_genome = true
) TO STDOUT CSV
