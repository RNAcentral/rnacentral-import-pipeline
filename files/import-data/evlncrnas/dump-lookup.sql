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
  -- Entries are matched by exact gene name and taxid, so every database that
  -- carries curated gene names for an EVlncRNAs species belongs here.
  AND upper(rnc_database.descr) IN (
    'ENSEMBL',
    'ENSEMBL_GENCODE',
    'ENSEMBL_PLANTS',
    'ENSEMBL_METAZOA',
    'ENSEMBL_FUNGI',
    'ENSEMBL_PROTISTS',
    'HGNC',
    'MGI',
    'RGD',
    'ZFIN',
    'FLYBASE',
    'WORMBASE',
    'SGD',
    'POMBASE',
    'TAIR',
    'REFSEQ',
    'LNCIPEDIA',
    'LNCBOOK',
    'NONCODE',
    'PLNCDB',
    'MALACARDS',
    'GENECARDS'
  )
  AND (
    gene <> ''
    OR external_id <> ''
    OR gene_synonym <> ''
    OR optional_id <> ''
  )

  ) TO STDOUT CSV HEADER
