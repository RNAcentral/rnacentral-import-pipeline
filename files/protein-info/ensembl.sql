SELECT
  gene.stable_id,
  gene.description,
  gene_xref.display_label,
  external_synonym.`synonym`
FROM gene
JOIN xref gene_xref ON gene_xref.xref_id = gene.display_xref_id
LEFT JOIN `external_synonym` ON external_synonym.xref_id = gene_xref.xref_id
WHERE
  gene.biotype = 'protein_coding'
;
