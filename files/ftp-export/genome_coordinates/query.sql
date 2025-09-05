COPY (
SELECT json_data
FROM (
  -- Transcript records
  SELECT
    json_build_object(
        'assembly_id', :'assembly_id',
        'region_id', regions.region_name,
        'rna_id', pre.id,
        'description', pre.short_description,
        'rna_type', pre.rna_type,
        'databases', regexp_split_to_array(pre."databases", ','),
        'providing_databases', COALESCE(array_agg(ac.database) FILTER (WHERE ac.database IS NOT NULL), '{}'),
        'chromosome', regions.chromosome,
        'strand', regions.strand,
        'identity', regions.identity,
        'was_mapped', regions.was_mapped,
        'exons', array_agg(distinct exons.*),
        'gene_name', genes.public_name,
        'gene_start', genes.start,
        'gene_end', genes.stop
    ) as json_data,
    regions.chromosome,
    regions.region_start,
    regions.id,
    'transcript'::text as record_type
  FROM rnc_rna_precomputed pre
  JOIN rnc_sequence_regions_active regions ON regions.urs_taxid = pre.id 
    AND regions.assembly_id = :'assembly_id'
  JOIN rnc_sequence_exons exons ON exons.region_id = regions.id
  LEFT JOIN rnc_gene_members gm ON gm.locus_id = regions.id
  LEFT JOIN rnc_genes genes ON genes.id = gm.rnc_gene_id
  LEFT JOIN rnc_accession_sequence_region sra ON sra.region_id = regions.id
  LEFT JOIN rnc_accessions ac ON sra.accession = ac.accession
  WHERE pre.is_active = true
  GROUP BY regions.id, regions.region_name, regions.chromosome, regions.strand, 
           regions.identity, regions.was_mapped, regions.region_start, pre.id, 
           pre.short_description, pre.rna_type, pre.databases, genes.public_name, genes.start, genes.stop

  UNION ALL

  -- Gene records
  SELECT
    json_build_object(
        'assembly_id', g.assembly_id,
        'region_id', g.public_name,
        'rna_id', NULL,
        'description', md.description,
        'rna_type', md.so_rna_type,
        'databases', ARRAY['RNAcentral-gene-prediction'],
        'providing_databases', NULL,
        'chromosome', g.chromosome,
        'strand', g.strand,
        'identity', NULL,
        'was_mapped', FALSE,
        'exons', NULL,
        'gene_name', g.public_name,
        'gene_start', g.start,
        'gene_end', g.stop
    ) as json_data,
    g.chromosome,
    g.start as region_start,
    g.id,
    'gene'::text as record_type
  FROM rnc_genes g
  JOIN rnc_gene_metadata md ON md.rnc_gene_id = g.id
  WHERE g.assembly_id = :'assembly_id'

  -- Sort everything together
  ORDER BY chromosome, region_start, record_type DESC, id
) combined_results
) TO STDOUT
