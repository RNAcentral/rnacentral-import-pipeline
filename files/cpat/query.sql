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
    -- Whole Ensembl mRNAs (dbid 61) are left out, their UTRs are kept. Read
    -- xref: precompute has not filled pre.databases for new sequences yet.
    AND NOT EXISTS (
      SELECT 1
      FROM xref
      JOIN rnc_accessions acc ON acc.accession = xref.ac
      WHERE xref.dbid = 61
        AND xref.deleted = 'N'
        AND xref.urs = pre.urs
        AND xref.taxid = pre.taxid
        AND acc.rna_type = 'SO:0000234'
    )
    AND pre.taxid = :taxid
    AND NOT exists(select 1 from cpat_results track where track.urs_taxid = pre.urs_taxid)
) TO STDOUT
