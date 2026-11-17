BEGIN;

-- xref is PARTITION BY LIST (dbid) and only prunes on a literal, so work out
-- which dbids this load touches (rnc_accessions.accession is its PK, so this is
-- an index probe per accession) and splice them into the join below. Same
-- pattern, and the full explanation, as 001__related-sequences.sql.
CREATE TEMP TABLE tmp_disease_dbids AS
SELECT DISTINCT d.id AS dbid
FROM (SELECT DISTINCT accession FROM load_rnc_disease_annotations) l
JOIN rnc_accessions acc ON acc.accession = l.accession
JOIN rnc_database d ON d.descr = acc.database;

-- Resolve accession -> urs_taxid. Annotations are loaded by accession because
-- the sequences may only have got a URS in this release.
DO $$
DECLARE
  v_unresolved bigint;
  v_unmatched_pmids bigint;
BEGIN
  EXECUTE format($q$
    CREATE TEMP TABLE tmp_disease_annotations AS
    SELECT DISTINCT
      xref.urs_taxid,
      load.accession,
      load.assigned_by,
      load.disease_text,
      nullif(load.ontology_term_id, '') AS ontology_term_id,
      load.match_type,
      nullif(load.pubmed_id, '') AS pubmed_id
    FROM load_rnc_disease_annotations load
    JOIN xref ON xref.ac = load.accession
    WHERE
      xref.deleted = 'N'
      AND xref.dbid = ANY(%L::int[])
  $q$, (SELECT coalesce(array_agg(dbid), '{}') FROM tmp_disease_dbids));

  -- Unresolved accessions are skipped, not an error: an entry can be parsed
  -- and then not make it into the release.
  SELECT count(DISTINCT load.accession) INTO v_unresolved
  FROM load_rnc_disease_annotations load
  WHERE NOT EXISTS (
    SELECT 1 FROM tmp_disease_annotations t WHERE t.accession = load.accession
  );
  RAISE NOTICE 'Skipped % load_rnc_disease_annotations accessions with no active xref', v_unresolved;

  SELECT count(DISTINCT t.pubmed_id) INTO v_unmatched_pmids
  FROM tmp_disease_annotations t
  WHERE
    t.pubmed_id IS NOT NULL
    AND NOT EXISTS (SELECT 1 FROM rnc_references refs WHERE refs.pmid = t.pubmed_id);
  RAISE NOTICE 'No rnc_references row for % disease annotation PMIDs, those links are not loaded', v_unmatched_pmids;
END $$;

-- Each load fully replaces its sources' annotations. Publication links go by
-- ON DELETE CASCADE. Keyed on the resolved rows, not the raw load: if nothing
-- for a source resolved to a URS, keep what it had rather than wipe it.
DELETE FROM rnc_disease_annotations
WHERE assigned_by IN (SELECT DISTINCT assigned_by FROM tmp_disease_annotations);

-- Deliberately no ON CONFLICT: one (urs_taxid, accession, disease_text) with
-- two different terms means the mapping is inconsistent, and that should fail
-- rather than pick one.
INSERT INTO rnc_disease_annotations (
    urs_taxid,
    accession,
    assigned_by,
    disease_text,
    ontology_term_id,
    match_type
) (
SELECT DISTINCT
    urs_taxid,
    accession,
    assigned_by,
    disease_text,
    ontology_term_id,
    match_type
FROM tmp_disease_annotations
);

INSERT INTO rnc_disease_annotation_publications (
    disease_annotation_id,
    reference_id
) (
SELECT DISTINCT
    annotations.disease_annotation_id,
    refs.id
FROM tmp_disease_annotations load
JOIN rnc_disease_annotations annotations
ON
    annotations.urs_taxid = load.urs_taxid
    AND annotations.accession = load.accession
    AND annotations.disease_text = load.disease_text
JOIN rnc_references refs
ON
    refs.pmid = load.pubmed_id
)
ON CONFLICT (disease_annotation_id, reference_id)
DO NOTHING
;

DROP TABLE load_rnc_disease_annotations;

COMMIT;
