LOAD CSV
FROM file 'data.csv'
HAVING FIELDS (database_name, task_name)
INTO {{PGDATABASE}}?load_ensembl_analysis_status
TARGET COLUMNS (database_name, task_name)
WITH fields terminated by ','

BEFORE LOAD DO
$$
drop table if exists load_ensembl_analysis_status;
$$,
$$
create table load_ensembl_analysis_status (
    database_name text,
    task_name text
);
$$

AFTER LOAD DO
$$
insert into ensembl_import_tracking (
    database_name text,
    task_name text,
    status
) (
select
    database_name text,
    task_name text,
    true
from load_ensembl_analysis_status
) ON CONFLICT DO UPDATE
set
    status = excluded.status
;
$$,
$$
drop table load_ensembl_analysis_status;
$$
;
