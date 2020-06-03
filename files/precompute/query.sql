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
        'species', coalesce(tax.name, acc.species),
        'common_name', coalesce(tax.common_name, acc.common_name),
        'feature_name', acc.feature_name,
        'ncrna_class', acc.ncrna_class,
        'locus_tag', acc.locus_tag,
        'organelle', acc.organelle,
        'lineage', coalesce(tax.lineage, acc.classification),
        'all_species', ARRAY[tax.name, acc.species::text],
        'all_common_names', ARRAY[tax.common_name, acc.common_name::text],
        'rna_type', acc.rna_type
      )
    ),
    'has_coordinates', exists(
        select 1
        from rnc_sequence_regions
        where urs_taxid = rna.upi || '_' || xref.taxid
    ),
    'deleted', array_agg(distinct xref.deleted = 'Y'),
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
    'last_release', max(xref.last),
    'chromosomes', ARRAY(
      select
        distinct chromosome
      from rnc_sequence_regions
        where urs_taxid = rna.upi || '_' || xref.taxid
    )
)
FROM rna
join :tablename todo on todo.upi = rna.upi
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
left join rnc_taxonomy tax on tax.id = xref.taxid
where
    todo.id BETWEEN :min AND :max
GROUP BY rna.upi, xref.taxid
ORDER BY rna.upi, xref.taxid
) TO STDOUT
