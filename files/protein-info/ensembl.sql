select
	gene.stable_id,
	gene.description,
	gene_xref.display_label,
	external_synonym.`synonym`

from gene
join xref gene_xref on gene_xref.xref_id = gene.display_xref_id
join `external_synonym` on external_synonym.xref_id = gene_xref.xref_id
where
	gene.biotype = 'protein_coding'
	and gene.stable_id = 'ENSG00000198901'
FIELDS TERMINATED BY ','
ENCLOSED BY '"'
LINES TERMINATED BY '\n'
;
