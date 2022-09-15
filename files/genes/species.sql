COPY (
  select
    distinct assembly_id, taxid
  from ensembl_assembly
  where selected_genome = true
) TO STDOUT CSV
