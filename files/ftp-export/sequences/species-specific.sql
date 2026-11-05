COPY (
-- Grouped by taxid so the by-species split can stream the result without sorting it.
SELECT DISTINCT ON (pre.taxid, pre.urs_taxid)
    json_build_object(
      'id', pre.urs_taxid,
      'description', pre.description,
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
FROM rnc_rna_precomputed pre
JOIN rna ON rna.urs = pre.urs
JOIN xref ON xref.urs = rna.urs AND xref.taxid = pre.taxid
WHERE
    pre.is_active = true
ORDER BY pre.taxid, pre.urs_taxid
) TO STDOUT
