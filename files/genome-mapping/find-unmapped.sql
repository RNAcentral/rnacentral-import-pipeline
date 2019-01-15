COPY (
    SELECT
    json_build_object(
      'id', pre.id,
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
    FROM rna
    JOIN rnc_rna_precomputed pre ON pre.upi = rna.upi
    LEFT JOIN rnc_sequence_regions regions
    ON
      regions.urs_taxid = pre.id
      AND regions.assembly_id = :'assembly_id'
      AND regions.was_mapped = false
    WHERE
      pre.taxid = :taxid
      AND regions.id IS NULL
      AND pre.is_active = true
) TO STDOUT
