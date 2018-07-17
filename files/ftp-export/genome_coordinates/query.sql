COPY (
  SELECT
    json_build_object(
      'rnacentral_id', concat_ws('_', t1.upi, t1.taxid)
      'chromosome', t2.name,
      'region_id', t2.accession,
      'strand', min(t2.strand),
      'exons', array_agg(
        json_build_object(
          'start', t2.primary_start,
          'stop', t2.primary_end
        )
      ),
      'rna_type', t3.rna_type,
      'databases', t3.databases
    )
  FROM xref t1
  JOIN rnc_coordinates t2
  ON
    t1.ac = t2.accession
  JOIN rnc_rna_precomputed t3
  ON
    t1.upi = t3.upi
    AND t1.taxid = t3.taxid
  WHERE t2.name IS NOT NULL
    AND t1.taxid = :taxid
    AND t1.deleted = 'N'
  GROUP BY accession, chromosome
  ORDER BY accession, primary_start, strand
UNION
  SELECT
    json_build_object(
      'rnacentral_id', t1.rna_id,
      'chromosome', chromosome,
      'region_id', region_id,
      'strand', min(t2.strand),
      'exons', array_agg(
        json_build_object(
          'start', start,
          'stop', stop
        )
      )
      'rna_type', t2.rna_type,
      'databases', t2.databases
    )
  FROM rnc_genome_mapping t1
  JOIN rnc_rna_precomputed t2
  ON t1.rna_id = t2.id
  WHERE
    t1.assembly_id = :assembly
  GROUP BY region_id
  ORDER BY region_id, start, strand
) TO STDOUT
