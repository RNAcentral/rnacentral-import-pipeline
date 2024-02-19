COPY(
  SELECT
    xref.upi as urs,
    xref.taxid as taxid,
    gene || '|' || external_id || '|' || gene_synonym || '|' || optional_id  as external_id
  FROM rnc_accessions
  JOIN xref
  ON xref.ac = rnc_accessions.accession
  WHERE xref.deleted = 'N'
  AND xref.dbid in (15,16,18,20,24,25,33,40,41,51)

  ) TO STDOUT CSV HEADER
