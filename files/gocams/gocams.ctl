LOAD CSV
FROM ALL FILENAMES MATCHING ~<gocams-memberships.*csv$>
WITH ENCODING UTF-8

HAVING FIELDS (
  urs_taxid,
  go_cam_id
)
INTO {{PGDATABASE}}?load_go_cam_memberships
TARGET COLUMNS (
  urs_taxid,
  go_cam_id
)

WITH
    fields escaped by double-quote,
    fields terminated by ','
;
