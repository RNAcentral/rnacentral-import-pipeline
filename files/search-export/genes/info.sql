COPY (
  SELECT
  json_build_object(
    'id', locus.id,
    'name', locus.locus_public_name,
    'assembly', locus.assembly
  )
  from rnc_locus locus
) TO STDOUT
