COPY (
  select
    distinct assembly_id, taxid
  from rnc_sequence_regions
) TO STDOUT
