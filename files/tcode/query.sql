COPY (
  SELECT
    json_build_object(
      'id', pre.urs_taxid,
      'sequence', coalesce(rna.seq_short, rna.seq_long)
    )
  FROM rnc_rna_precomputed pre
  JOIN rna
    ON rna.urs = pre.urs
  WHERE
    rna.id >= :min
    AND rna.id < :max
    AND pre.is_active = true
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
    AND pre.taxid IS NOT NULL
    AND NOT EXISTS (select 1 from tcode_results as tc where tc.urs_taxid = pre.urs_taxid)
    -- Filter sequences with non-ACGTUN letters or >4 Ns to avoid tcode crashes
    AND coalesce(rna.seq_short, rna.seq_long) ~ '^[ACGTUNacgtun]+$'
    AND coalesce(rna.seq_short, rna.seq_long) !~ '([Nn][^Nn]*){5}'
    AND rna.len > :min_len
) TO STDOUT
