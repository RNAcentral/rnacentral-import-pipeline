LOAD CSV
FROM 'imported.csv'
HAVING FIELDS (task_name, database_name)
INTO {{PGDATABASE}}?load_ensembl_analysis_status
TARGET COLUMNS (task_name, database_name)
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
    database_name,
    task_name,
    was_imported
) (
select
    database_name,
    task_name,
    true
from load_ensembl_analysis_status
) ON CONFLICT (database_name, task_name) DO UPDATE
set
    was_imported = excluded.was_imported
;
$$,
$$
drop table load_ensembl_analysis_status;
$$
;
