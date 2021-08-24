COPY (
  SELECT
  json_build_object(
    'id', locus.id,
    'name', locus.locus_public_name,
  )
  from rnc_locus locus
) TO STDOUT
