COPY (
  SELECT
    json_build_object(
      'assembly_id', genome.assembly_id,
      'assembly_full_name', genome.assembly_full_name,
      'gca_accession', genome.gca_accession,
      'assembly_ucsc', genome.assembly_ucsc,
      'common_name', genome.common_name,
      'taxid', genome.taxid,
      'ensembl_url', genome.ensembl_url,
      'division', genome.division,
      'blat_mapping', genome.blat_mapping,
      'example', json_build_object(
        'chromosome', genome.chromosome,
        'start', genome.assembly.start,
        'end', genome.assembly_end
      )
    )
    FROM ensembl_assembly genome
) TO STDOUT
