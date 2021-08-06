COPY (
  select
    distinct assembly_id, taxid
  from from ensembl_assembly
) TO STDOUT CSV
