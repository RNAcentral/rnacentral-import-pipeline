LOAD CSV
FROM file 'imported-data.csv'
HAVING FIELDS (task_name)
INTO {{PGDATABASE}}?load_ensembl_analysis_status
TARGET COLUMNS (task_name)
WITH fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_metadata_import_status;
$$,
$$
create table load_metadata_import_status (
    task_name text
);
$$

AFTER LOAD DO
$$
insert into metadata_import_status (
    task_name,
    status
) (
select
    task_name,
    true
from load_ensembl_analysis_status
) ON CONFLICT DO UPDATE
set
    status = excluded.status
;
$$,
$$
drop table load_metadata_import_status;
$$
;
