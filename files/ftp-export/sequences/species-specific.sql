COPY (
SELECT DISTINCT ON (pre.id)
    json_build_object(
      'id', pre.id,
      'description', pre.description,
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
FROM rnc_rna_precomputed pre
JOIN rna ON rna.upi = pre.upi
JOIN xref ON xref.upi = rna.upi AND xref.taxid = pre.taxid
WHERE
    pre.is_active = true
) TO STDOUT
