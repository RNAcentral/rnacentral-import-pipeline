COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'cross_references', array_agg(
          json_build_object(
            'name', acc."database",
            'external_id', acc.external_id,
            'optional_id', acc.optional_id,
            'accession', acc.accession,
            'non_coding_id', acc.non_coding_id,
            'parent_accession', acc.parent_ac || '.' || acc.seq_version
          )
      )
  )
  FROM :tablename todo
  JOIN xref xref
  ON
    xref.upi = todo.urs
  JOIN rnc_accessions acc
  ON
    xref.ac = acc.accession
  GROUP BY todo.id
) TO STDOUT
