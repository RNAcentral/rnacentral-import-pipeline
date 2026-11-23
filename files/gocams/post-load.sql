BEGIN;

-- One row per (RNA, GO-CAM): an RNA can sit in several models.
CREATE TABLE IF NOT EXISTS go_cam_memberships (
    urs_taxid varchar(50) NOT NULL,
    go_cam_id varchar(50) NOT NULL,
    PRIMARY KEY (urs_taxid, go_cam_id)
);

-- The primary key serves lookups by urs_taxid. This serves "which RNAs are in
-- this model" lookups.
CREATE INDEX IF NOT EXISTS go_cam_memberships_go_cam_id_idx
    ON go_cam_memberships (go_cam_id);

-- GO-CAM has no release cycle, so files/qc/analyze.sql reads this to see the
-- load ran.
CREATE TABLE IF NOT EXISTS pipeline_tracking_go_cam_memberships (
    source text PRIMARY KEY,
    last_run timestamptz NOT NULL,
    rows_loaded bigint NOT NULL
);

-- Full snapshot, not an upsert: an updated GO-CAM gets a new id and models can
-- be retracted, so old (urs_taxid, go_cam_id) pairs must go. TRUNCATE is
-- transactional: readers see the old rows until COMMIT, a failure rolls back.
DO $$
BEGIN
    IF NOT EXISTS (SELECT 1 FROM load_go_cam_memberships) THEN
        RAISE EXCEPTION 'load_go_cam_memberships is empty; refusing to wipe go_cam_memberships';
    END IF;
END
$$;

TRUNCATE go_cam_memberships;

INSERT INTO go_cam_memberships (urs_taxid, go_cam_id)
SELECT DISTINCT urs_taxid, go_cam_id
FROM load_go_cam_memberships
;

INSERT INTO pipeline_tracking_go_cam_memberships (source, last_run, rows_loaded)
VALUES ('go-cam', now(), (SELECT count(*) FROM load_go_cam_memberships))
ON CONFLICT (source) DO UPDATE
SET
    last_run = excluded.last_run,
    rows_loaded = excluded.rows_loaded
;

DROP TABLE load_go_cam_memberships;

COMMIT;
