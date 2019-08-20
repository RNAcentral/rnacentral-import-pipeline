COPY (
SELECT
  json_build_object(
    'rna_id', pre.id,
    'secondary', json_build_object(
      'has_secondary', true,
      'secondary_structure_model', 'gtrnadb'
    )
  )
FROM rnc_rna_precomputed pre
join xref on xref.upi = pre.upi
join rnc_secondary_structure ss
on
  ss.rnc_accession_id = xref.ac
where
  pre.is_active = true
  and xref.dbid = 8
  and xref.deleted = 'N'
) TO STDOUT
