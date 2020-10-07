COPY (
  select
    distinct assembly_id
  from rnc_sequence_regions
)
