COPY (
select distinct
  assembly.assembly_id,
  assembly.ensembl_url,
  assembly.taxid
from ensembl_assembly assembly
join ensembl_coordinate_systems coord
ON
  coord.assembly_id = assembly.assembly_id
  and coord.is_reference = true
) TO STDOUT CSV
