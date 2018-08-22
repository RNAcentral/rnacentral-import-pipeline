COPY (
SELECT
    json_build_object(
      'id', pre.id,
      'description', pre.description,
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
FROM rnc_rna_precomputed pre
JOIN rna on rna.upi = pre.upi
WHERE
    pre.is_active = false
    AND pre.taxid is null
) TO STDOUT
