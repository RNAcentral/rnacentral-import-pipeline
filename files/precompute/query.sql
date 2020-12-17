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
        'so_rna_type', acc.rna_type
      )
    ),
    'coordinates', array_agg(json_build_object(
        'assembly_id', region.assembly_id,
        'chromosome', region.chromosome,
        'strand', region.strand,
        'start', region.region_start,
        'stop', region.region_stop
      )),
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
    'last_release', max(xref.last)
)
FROM rna
JOIN :tablename todo ON todo.upi = rna.upi
JOIN xref 
ON 
    xref.upi = rna.upi
JOIN rnc_accessions acc
ON
    acc.accession = xref.ac
LEFT JOIN rnc_rna_precomputed prev
ON
    prev.upi = rna.upi
    AND prev.taxid = xref.taxid
LEFT JOIN rnc_sequence_regions region 
ON 
    region.urs_taxid = xref.upi || '_' || xref.taxid
LEFT JOIN rfam_model_hits hits 
ON 
    hits.upi = xref.upi
LEFT JOIN rfam_models models 
ON 
    models.rfam_model_id = hits.rfam_model_id
JOIN rnc_database db
ON
    db.id = xref.dbid
LEFT JOIN rnc_taxonomy tax ON tax.id = xref.taxid
WHERE
    todo.id BETWEEN :min AND :max
GROUP BY rna.upi, xref.taxid
ORDER BY rna.upi, xref.taxid
) TO STDOUT
