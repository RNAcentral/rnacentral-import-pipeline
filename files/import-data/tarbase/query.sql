COPY (
SELECT DISTINCT ON (acc.optional_id)
  json_build_object(
    'id', acc.optional_id,
    'sequence', coalesce(rna.seq_short, rna.seq_long)
  )
FROM rnc_accessions acc
JOIN xref_p4_not_deleted xref ON xref.ac = acc.accession
JOIN rna on rna.upi = xref.upi
WHERE
  xref.deleted = 'N'
  AND xref.dbid = 4
  AND acc.database = 'MIRBASE'
) TO STDOUT
