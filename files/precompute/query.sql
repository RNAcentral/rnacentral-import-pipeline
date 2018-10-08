COPY (
select
  json_build_object(
    'upi', rna.upi,
    'taxid', xref.taxid,
    'length', rna.len,
    'accessions', array_agg(json_build_object(
        'is_active', xref.deleted = 'N',
        'description', acc.description,
        'gene', acc.gene,
        'optional_id', acc.optional_id,
        'pretty_database', db.display_name,
        'species', acc.species,
        'common_name', acc.common_name,
        'feature_name', acc.feature_name,
        'ncrna_class', acc.ncrna_class,
        'locus_tag', acc.locus_tag,
        'organelle', acc.organelle,
        'lineage', acc.classification
      )
    ),
    'deleted', array_agg(distinct xref.deleted = 'Y'),
    'xref_has_coordinates', array_agg(exists(
        select 1
        from rnc_coordinates coord
        where coord.accession = xref.ac
    )),
  'rna_was_mapped', exists(
    select 1
    from rnc_genome_mapping mapping
    where mapping.upi = rna.upi and mapping.taxid = xref.taxid
  ),
  'previous', array_agg(row_to_json(prev.*)),
  'hits', array_agg(json_build_object(
      'rfam_hit_id', hits.rfam_hit_id,
      'model', hits.rfam_model_id,
      'model_rna_type', models.rna_type,
      'model_domain', models.domain,
      'model_name', models.short_name,
      'model_long_name', models.long_name,
      'model_completeness', hits.model_completeness,
      'model_start', hits.model_start,
      'model_stop', hits.model_stop,
      'sequence_completeness', hits.sequence_completeness,
      'sequence_start', hits.sequence_start,
      'sequence_stop', hits.sequence_stop
  )),
  'last_release', max(xref.last)
)
FROM rna
join upis_to_precompute todo on todo.upi = rna.upi
join xref on xref.upi = rna.upi
join rnc_accessions acc
on
    acc.accession = xref.ac
left join rnc_rna_precomputed prev
on
    prev.upi = rna.upi
    and prev.taxid = xref.taxid
left join rfam_model_hits hits ON hits.upi = xref.upi
left join rfam_models models ON models.rfam_model_id = hits.rfam_model_id
join rnc_database db
ON
    db.id = xref.dbid
where
    todo.id BETWEEN :min AND :max
GROUP BY rna.upi, xref.taxid
ORDER BY rna.upi, xref.taxid
) TO STDOUT
