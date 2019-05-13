COPY (
SELECT
  optional_id
FROM rnc_accessions acc
JOIN xref on xref.ac = acc.accession
WHERE
    xref.deleted = 'N'
    AND xref.dbid = 8
) TO STDOUT
