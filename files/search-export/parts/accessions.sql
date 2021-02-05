COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid
      'accession', todo.accession,
      'common_name', todo.common_name,
      'database', todo.database,
      'external_id', todo.external_id,
      'functions', todo.function,
      'gene_synonyms', todo.gene_synonym,
      'genes', todo.gene,
      'locus_tags', todo.locus_tag,
      'non_coding_id', todo.non_coding_id,
      'notes', todo.note,
      'optional_id', todo.optional_id,
      'organelles', todo.organelle,
      'parent_accession', todo.parent_acccession
      'products', todo.product,
      'species', todo.species,
      'standard_names', todo.standard_name,
      'tax_strings', todo.classification,
    )
  FROM search_export_accessions todo
  ORDER BY todo.search_export_id
) TO STDOUT
