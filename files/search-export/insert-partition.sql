INSERT INTO search_export_accessions (
  search_export_id,
  urs_taxid,
  accession,
  classification,
  common_name,
  database,
  external_id,
  function,
  gene,
  gene_synonym,
  genes,
  locus_tag,
  locus_tags,
  non_coding_id,
  note,
  optional_id,
  organelle,
  organelles,
  parent_accession,
  product,
  species,
  standard_name
) (
SELECT
  todo.id,
  todo.urs_taxid,
  acc.accession,
  acc.classification,
  acc.common_name,
  acc.database,
  acc.external_id,
  acc.function,
  acc.gene,
  acc.gene_synonym,
  acc.genes,
  acc.locus_tag,
  acc.locus_tags,
  acc.non_coding_id,
  acc.note,
  acc.optional_id,
  acc.organelle,
  acc.organelles,
  acc.parent_ac || '.' || acc.seq_version,
  acc.product,
  acc.species,
  acc.standard_name
FROM search_export_urs todo
JOIN :partition xref
ON
  xref.urs = todo.urs 
  AND xref.taxid = todo.taxid
JOIN rnc_accessions acc
ON
  acc.accession = xref.ac
);
