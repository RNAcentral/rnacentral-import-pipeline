SELECT meta_key, meta_value
FROM meta
WHERE meta_key IN (
  'assembly.accession',
  'assembly.default',
  'assembly.long_name',
  'assembly.name',
  'assembly.ucsc_alias',
  'species.division',
  'species.production_name',
  'species.taxonomy_id',
  'species.common_name',
  'species.scientific_name',
    'species.url',
    'species.division'
)

