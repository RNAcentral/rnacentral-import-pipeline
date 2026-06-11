COPY (
SELECT
    json_build_object(
      'id', pre.upi,
      'description',
        (CASE WHEN count(DISTINCT pre.rna_type) = 1 THEN max(pre.rna_type) ELSE 'ncRNA' END)
        || ' found in ' || count(DISTINCT pre.taxid) || ' taxa',
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
FROM rnc_rna_precomputed pre
JOIN rna on rna.upi = pre.upi
GROUP BY pre.upi, rna.seq_short, rna.seq_long
HAVING bool_or(pre.is_active) = false
) TO STDOUT
