COPY (
  SELECT
    json_build_object(
      'id', xref.upi || '_' || xref.taxid,
      'name', acc."database",
      'external_id', acc.external_id,
      'optional_id', acc.optional_id,
      'accession', acc.accession,
      'non_coding_id', acc.non_coding_id,
      'parent_accession', acc.parent_ac || '.' || acc.seq_version
  )
  FROM xref xref
  JOIN rnc_accessions acc
  ON
    xref.ac = acc.accession
) TO STDOUT
