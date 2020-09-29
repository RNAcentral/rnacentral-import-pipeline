COPY (
  SELECT
    json_build_object(
      'id', xref.upi || '_' || xref.taxid,
      'species', acc.species,
      'organelles', acc.organelle,
      'product', acc.product,
      'tax_strings', acc.classification,
      'functions', acc.function,
      'genes', acc.gene,
      'gene_synonyms', acc.gene_synonym,
      'common_name', acc.common_name,
      'notes', acc.note,
      'locus_tags', acc.locus_tag,
      'standard_names', acc.standard_name,
      'products', acc.product
    )
  FROM xref
  ON
    xref.upi = todo.urs
  JOIN rnc_accessions acc
  ON
    xref.ac = acc.accession
) TO STDOUT
