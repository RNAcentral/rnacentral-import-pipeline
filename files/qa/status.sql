COPY (
  select
    json_build_object(
      'upi', xref.upi,
      'taxid', xref.taxid,
      'hits', array_agg(json_build_object(
          'model', hits.rfam_model_id,
          'model_rna_type', models.rna_type,
          'model_domain', models.domain,
          'model_completeness', hits.model_completeness,
          'sequence_completeness', hits.sequence_completeness
        )),
      'rna_types', array_agg(distinct acc.
      'organelles', array_agg(distinct acc.organelle),
      'descriptions', array_agg(distinct acc.description),
      'lineages', array_agg(distinct acc.classification)
    )
  from xref
  join rnc_accessions acc on acc.accession = xref.accession
  left join rfam_model_hits hits on hits.upi = xref.upi
  left join rfam_models models on models.rfam_model_id = hits.rfam_model_id
  group by xref.upi, xref.taxid
) TO STDOUT
