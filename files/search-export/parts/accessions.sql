COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'species', array_agg(acc.species),
      'organelles', array_agg(acc.organelle),
      'product', array_agg(acc.product),
      'tax_strings', array_agg(acc.classification),
      'functions', array_agg(acc.function),
      'genes', array_agg(acc.gene),
      'gene_synonyms', array_agg(acc.gene_synonym),
      'common_name', array_agg(acc.common_name),
      'notes', array_agg(acc.note),
      'locus_tags', array_agg(acc.locus_tag),
      'standard_names', array_agg(acc.standard_name),
      'products', array_agg(acc.product)
    )
  FROM :tablename todo
  JOIN xref xref
  ON
    xref.upi = todo.urs
  JOIN rnc_accessions acc
  ON
    xref.ac = acc.accession
  GROUP BY todo.id
) TO STDOUT
