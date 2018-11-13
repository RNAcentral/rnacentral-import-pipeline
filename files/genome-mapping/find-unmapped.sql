COPY (
    select
    json_build_object(
      'id', pre.id,
      'sequence', COALESCE(rna.seq_short, rna.seq_long)
    )
    FROM rna
    JOIN rnc_rna_precomputed pre on pre.upi = rna.upi
    LEFT JOIN rnc_sequence_regions regions
    ON
      regions.urs_taxid = pre.id
      and regions.assembly_id = :'assembly_id'
    WHERE
      pre.taxid = :taxid
      and regions.id is null
      and pre.is_active = true
) TO STDOUT
