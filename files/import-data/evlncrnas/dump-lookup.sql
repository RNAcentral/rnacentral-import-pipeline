COPY(
  SELECT
    xref.urs as urs,
    xref.taxid as taxid,
    gene || '|' || external_id || '|' || gene_synonym || '|' || optional_id  as external_id
  FROM rnc_accessions
  JOIN xref
  ON xref.ac = rnc_accessions.accession
  JOIN rnc_database
  ON rnc_database.id = xref.dbid
  WHERE xref.deleted = 'N'
  AND upper(rnc_database.descr) IN (
    'WORMBASE',
    'SGD',
    'POMBASE',
    'LNCIPEDIA',
    'FLYBASE',
    'ENSEMBL',
    'LNCBOOK',
    'MALACARDS',
    'GENECARDS',
    'EXPRESSION_ATLAS'
  )

  ) TO STDOUT CSV HEADER
