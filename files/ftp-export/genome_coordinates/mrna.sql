-- Ensembl protein coding mRNAs on one assembly, for the IGV gene-model track.
-- Exported as CSV so the JSON survives the note's quotes.
COPY (
SELECT
  json_build_object(
    'transcript', acc.accession,
    'rna_id', xref.urs || '_' || xref.taxid,
    'gene', acc.optional_id,
    'gene_name', acc.locus_tag,
    'note', acc.note,
    -- A transcript's UTRs are imported as <accession>:<five|three>_prime_UTR
    'five_prime_utr', (
      SELECT utr.urs || '_' || utr.taxid FROM xref utr
      WHERE utr.dbid = 61 AND utr.deleted = 'N'
        AND utr.ac = acc.accession || ':five_prime_UTR'
    ),
    'three_prime_utr', (
      SELECT utr.urs || '_' || utr.taxid FROM xref utr
      WHERE utr.dbid = 61 AND utr.deleted = 'N'
        AND utr.ac = acc.accession || ':three_prime_UTR'
    ),
    'chromosome', regions.chromosome,
    'strand', regions.strand,
    'exons', array_agg(
      json_build_object('start', exons.exon_start, 'stop', exons.exon_stop)
      ORDER BY exons.exon_start
    )
  )
FROM xref
JOIN rnc_accessions acc ON acc.accession = xref.ac
JOIN rnc_accession_sequence_region sra ON sra.accession = acc.accession
JOIN rnc_sequence_regions regions ON regions.id = sra.region_id
JOIN rnc_sequence_exons exons ON exons.region_id = regions.id
WHERE
  xref.dbid = 61
  AND xref.deleted = 'N'
  AND acc.feature_name = 'mRNA'
  AND regions.assembly_id = :'assembly_id'
GROUP BY
  acc.accession, xref.urs, xref.taxid, acc.optional_id, acc.locus_tag,
  acc.note, regions.chromosome, regions.strand
) TO STDOUT CSV
