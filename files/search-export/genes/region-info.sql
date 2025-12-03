COPY (
  SELECT
  json_build_object(
    'id', todo.id,
    'urs_taxid', todo.urs_taxid,
    'gene_id', gene.id,
    'gene_name', gene.public_name,
    'gene_description', metadata.description,
    'so_rna_type', metadata.so_rna_type,
    'short_description', metadata.short_description,
    'start', gene."start",
    'stop', gene.stop,
    'assembly_id', gene.assembly_id,
    'member_count', gene.member_count
  )
  FROM search_export_urs todo
  JOIN rnc_sequence_regions regions
  ON
    regions.urs_taxid = todo.urs_taxid
  JOIN rnc_gene_members mem
  ON
    mem.locus_id = regions.id
  JOIN rnc_genes gene
  ON
    gene.id = mem.rnc_gene_id
  JOIN rnc_gene_metadata metadata
  ON
    metadata.rnc_gene_id = gene.id
  ORDER BY todo.id
) TO STDOUT
