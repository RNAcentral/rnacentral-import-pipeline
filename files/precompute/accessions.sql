COPY (
select
  json_build_object(
    'id', todo.id,
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
    )
FROM :tablename todo
) TO STDOUT
