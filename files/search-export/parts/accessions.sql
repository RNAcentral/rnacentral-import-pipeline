COPY (
  SELECT
    json_build_object(
      'urs_taxid', todo.urs_taxid
      'species', todo.species,
      'organelles', todo.organelle,
      'product', todo.product,
      'tax_strings', todo.classification,
      'functions', todo.function,
      'genes', todo.gene,
      'gene_synonyms', todo.gene_synonym,
      'common_name', todo.common_name,
      'notes', todo.note,
      'locus_tags', todo.locus_tag,
      'standard_names', todo.standard_name,
      'products', todo.product,
      'database', todo.database,
      'external_id', todo.external_id,
      'optional_id', todo.optional_id,
      'accession', todo.accession,
      'non_coding_id', todo.non_coding_id,
      'parent_accession', todo.parent_acccession
    )
  FROM search_export_accessions todo
  ORDER BY todo.search_export_id
) TO STDOUT
