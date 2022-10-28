COPY(
  SELECT urs_taxid,
    xref.taxid as taxid,
    gene || '|' || external_id || '|' || gene_synonym || '|' || optional_id  as external_id,
    description,
    seq_version,
    assembly_id,
    region_start,
    region_stop,
    rsr.chromosome,
    strand,
    rna_type,
    COALESCE(seq_short, seq_long) as seq
  FROM rnc_accessions
  JOIN xref
  ON xref.ac = rnc_accessions.accession

  JOIN rna
  ON xref.upi = rna.upi

  JOIN rnc_accession_sequence_region rasr
  ON rasr.accession = xref.ac

  JOIN rnc_sequence_regions rsr
  ON rsr.id = region_id

  ) TO STDOUT CSV HEADER
