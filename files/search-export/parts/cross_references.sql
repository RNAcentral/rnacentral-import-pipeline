COPY (
  SELECT
    json_build_object(
      'id', xref.upi || '_' || xref.taxid,
      'accession', acc.accession,
      'common_name', acc.common_name,
      'database', acc."database",
      'external_id', acc.external_id,
      'functions', acc.function,
      'gene_synonyms', acc.gene_synonym,
      'genes', acc.gene,
      'locus_tags', acc.locus_tag,
      'non_coding_id', acc.non_coding_id,
      'notes', acc.note,
      'optional_id', acc.optional_id,
      'organelles', acc.organelle,
      'parent_accession', acc.parent_ac || '.' || acc.seq_version,
      'product', acc.product,
      'species', acc.species,
      'standard_names', acc.standard_name,
      'tax_strings', acc.classification
  )
  FROM xref xref
  JOIN rnc_accessions acc
  ON
    xref.ac = acc.accession
  WHERE
    xref.deleted = 'N'
) TO STDOUT
