COPY (
SELECT
  json_build_object(
    'taxid', taxid,
    'model_name', model_name,
    'rna_type', rna_type,
    'so_term_id', so_term_id
  )
from rnc_secondary_structure_layout_models
where
  model_source = 'crw'
) TO STDOUT
