COPY (
SELECT
    json_build_object(
      'id', pre.urs,
      'description',
        (CASE WHEN count(DISTINCT pre.rna_type) = 1 THEN max(pre.rna_type) ELSE 'ncRNA' END)
        || ' found in ' || count(DISTINCT pre.taxid) || ' taxa',
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
FROM rnc_rna_precomputed pre
JOIN rna ON rna.urs = pre.urs
WHERE
    pre.is_active = true
GROUP BY pre.urs, rna.seq_short, rna.seq_long
) TO STDOUT
