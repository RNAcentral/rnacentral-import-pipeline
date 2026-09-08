COPY (
SELECT
    json_build_object(
      'id', pre.urs_taxid,
      'description', pre.description,
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
FROM rnc_rna_precomputed pre
JOIN rna ON rna.urs = pre.urs
WHERE
    pre.is_active = true
    AND pre.taxid is not null
    AND pre.databases ilike :'db'
ORDER BY pre.urs_taxid
) TO STDOUT
