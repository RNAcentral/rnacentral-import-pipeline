-- Run this as a one-off the first time we load the gocams

BEGIN;

CREATE TABLE IF NOT EXISTS rnacen.go_cam_memberships (
    urs_taxid varchar(50) NOT NULL,
    go_cam_id varchar(50) NOT NULL,
    PRIMARY KEY (urs_taxid, go_cam_id)
);

-- The primary key serves lookups by urs_taxid. This serves "which RNAs are in
-- this model" lookups.
CREATE INDEX IF NOT EXISTS go_cam_memberships_go_cam_id_idx
    ON rnacen.go_cam_memberships (go_cam_id);

CREATE TABLE IF NOT EXISTS rnacen.pipeline_tracking_go_cam_memberships (
    source text PRIMARY KEY,
    last_run timestamptz NOT NULL,
    rows_loaded bigint NOT NULL
);

COMMIT;
