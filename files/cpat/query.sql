COPY (
  SELECT
    json_build_object(
      'id', pre.urs_taxid,
      'sequence', coalesce(rna.seq_short, rna.seq_long)
    )
  FROM rnc_rna_precomputed pre
  join rna on rna.urs = pre.urs
  where
    pre.is_active = true
    AND pre.taxid = :taxid
    AND NOT exists(select 1 from cpat_results track where track.urs_taxid = pre.urs_taxid)
) TO STDOUT
