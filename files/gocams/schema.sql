DROP TABLE IF EXISTS load_go_cam_memberships;
CREATE UNLOGGED TABLE load_go_cam_memberships (
    urs_taxid varchar(50),
    go_cam_id varchar(50)
);
