COPY (
SELECT
  json_build_object(
      'group', case
      when acc.database = 'ENSEMBL' then acc.gene,
      end,
      'urs_taxid', max(pre.id),
      'rna_type',  max(pre.rna_type),
      'exons', array_agg(distinct exons.*)
  )
FROM rnc_rna_precomputed pre
JOIN rnc_sequence_regions regions
ON
  regions.urs_taxid = pre.id
JOIN rnc_sequence_exons exons
ON
  exons.region_id = regions.id
JOIN xref 
ON
  xref.upi = pre.upi
  and xref.taxid = pre.taxid
JOIN rnc_accessions acc
ON
  acc.accession = xref.ac
WHERE
  pre.is_active = true
  AND (
    (acc.database = 'ENSEMBL' 
      AND acc.gene IN (
        'ENSG00000229807.13',
        'ENSG00000270641.1',
      )
  )
  )
GROUP BY regions.id
ORDER BY max(regions.chromosome), regions.region_start, regions.id
) TO STDOUT


-- Get coordaintes for:
-- /URS0000735371/9606 all in that region
-- URS00004B3771/9606 all in that region
-- ENSG00000229807.13 all XIST transcripts
-- URS000075E149/9606 everything overlapping
-- URS0000ABD82A/9606
-- FBgn0267665
-- FBgn0287460
