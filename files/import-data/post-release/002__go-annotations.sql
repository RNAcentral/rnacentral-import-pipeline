BEGIN;

-- QuickGO has no release cycle (no rnc_release/release_stats row), so this is
-- how the import QC report finds out it ran at all.
CREATE TABLE IF NOT EXISTS pipeline_tracking_go_annotations (
    source text PRIMARY KEY,
    last_run timestamptz NOT NULL,
    rows_loaded bigint NOT NULL
);

INSERT INTO go_term_annotations (
    urs_taxid,
    qualifier,
    ontology_term_id,
    evidence_code,
    assigned_by,
    extensions
) (
SELECT
    rna_id,
    qualifier,
    ontology_term_id,
    evidence_code,
    assigned_by,
    extensions
FROM load_go_term_annotations
)
ON CONFLICT (urs_taxid, qualifier, ontology_term_id, evidence_code, assigned_by)
DO UPDATE
SET
    extensions = excluded.extensions
;

INSERT INTO pipeline_tracking_go_annotations (source, last_run, rows_loaded)
VALUES ('quickgo', now(), (SELECT count(*) FROM load_go_term_annotations))
ON CONFLICT (source) DO UPDATE
SET
    last_run = excluded.last_run,
    rows_loaded = excluded.rows_loaded
;

DROP TABLE load_go_term_annotations;

COMMIT;
