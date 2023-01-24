COPY(
  SELECT
    xref.upi as urs,
    xref.taxid as taxid,
    gene || '|' || external_id || '|' || gene_synonym || '|' || optional_id  as external_id
  FROM rnc_accessions
  JOIN xref
  ON xref.ac = rnc_accessions.accession
  WHERE xref.deleted = 'N'
  AND xref.dbid in (15,16,18,20,24,25,33,40,41,51)

  ) TO STDOUT CSV HEADER

-- COPY(
--   SELECT urs_taxid,
--     xref.taxid as taxid,
--     gene || '|' || external_id || '|' || gene_synonym || '|' || optional_id  as external_id,
--     description,
--     seq_version,
--     assembly_id,
--     region_start,
--     region_stop,
--     rsr.chromosome,
--     strand,
--     rna_type,
--     COALESCE(seq_short, seq_long) as seq
--   FROM rnc_accessions
--   JOIN xref
--   ON xref.ac = rnc_accessions.accession

--   JOIN rna
--   ON xref.upi = rna.upi

--   JOIN rnc_accession_sequence_region rasr
--   ON rasr.accession = xref.ac

--   JOIN rnc_sequence_regions rsr
--   ON rsr.id = region_id

--   WHERE xref.deleted = 'N'
--   ) TO STDOUT CSV HEADER

-- COPY(
-- select
-- 	xref.upi || '_' || xref.taxid as urs_taxid,
-- 	ra.optional_id,
-- 	external_id,
--     exon_start,
--     exon_stop,
--     strand,
--     sr.chromosome
-- from xref join
-- rnc_accessions ra on
-- ra.accession = xref.ac
-- join rnc_sequence_regions sr
-- on xref.upi || '_' || xref.taxid = sr.urs_taxid
-- join rnc_sequence_exons ex on ex.region_id = sr.id
-- where xref.deleted = 'N'
--  ) TO STDOUT CSV HEADER
