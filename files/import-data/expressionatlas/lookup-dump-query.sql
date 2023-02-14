COPY(
  SELECT xref.upi || '_' || xref.taxid as urs_taxid,
    xref.taxid as taxid,
    gene || '|' || external_id || '|' || gene_synonym || '|' || optional_id  as external_id,
    description,
    seq_version,
    rna_type,
    COALESCE(seq_short, seq_long) as seq
  FROM rnc_accessions
  JOIN xref
  ON xref.ac = rnc_accessions.accession

  JOIN rna
  ON xref.upi = rna.upi


  ) TO STDOUT CSV HEADER
