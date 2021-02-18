COPY (
  SELECT
    json_build_object(
      'id', todo.search_export_id,
      'urs_taxid', todo.urs_taxid,
      'accession', todo.accession,
      'common_name', COALESCE(tax.common_name, todo.common_name),
      'database', todo.database,
      'external_id', todo.external_id,
      'function', todo.function,
      'gene_synonyms', todo.gene_synonym,
      'gene', todo.gene,
      'locus_tag', todo.locus_tag,
      'non_coding_id', todo.non_coding_id,
      'notes', todo.note,
      'optional_id', todo.optional_id,
      'organelle', todo.organelle,
      'parent_accession', todo.parent_accession,
      'product', todo.product,
      'species', COALESCE(tax.name, todo.species),
      'standard_name', todo.standard_name,
      'tax_string', COALESCE(tax.lineage, todo.lineage),
      'authors', refs.authors,
      'journal', refs.location,
      'pub_title', refs.title,
      'pub_id', refs.id,
      'pubmed_id', refs.pmid,
      'doi', refs.doi
    )
  FROM search_export_accessions todo
  LEFT JOIN rnc_taxonomy tax on tax.id = todo.taxid
  LEFT JOIN rnc_reference_map ref_map
  ON
    ref_map.accession = todo.accession
  LEFT JOIN rnc_references refs
  ON
    refs.id = ref_map.reference_id
  WHERE
    todo.search_export_id BETWEEN :min AND :max
  ORDER BY todo.search_export_id
) TO STDOUT
