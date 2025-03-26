CREATE TEMP TABLE taxids_to_fetch (
  taxid bigint PRIMARY KEY
);

\copy taxids_to_fetch FROM 'taxids_to_fetch';


COPY(
  SELECT xref.upi || '_' || xref.taxid as urs_taxid,
    xref.taxid as taxid,
    split_part(gene, '.', 1) as gene,
    external_id,
    gene_synonym ,
    optional_id,
    description,
    seq_version,
    rna_type,
    COALESCE(seq_short, seq_long) as seq
  FROM rna
  JOIN xref
  ON xref.upi = rna.upi
  join rnc_accessions
  ON xref.ac = rnc_accessions.accession

  WHERE xref.deleted = 'N'
  AND xref.dbid IN (25, 31, 34, 35, 36, 47)
  AND (
    gene <> ''
    OR
    external_id <> ''
    OR
    gene_synonym <> ''
    OR optional_id <> ''
  )
  AND xref.taxid IN (SELECT taxid FROM taxids_to_fetch)


  ) TO STDOUT CSV HEADER
