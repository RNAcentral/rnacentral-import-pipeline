COPY (
  select
    distinct assembly_id, taxid
  from ensembl_assembly
) TO STDOUT CSV
