LOAD CSV
FROM ALL FILENAMES MATCHING ~<ensembl_imported.*csv$>
HAVING FIELDS (task_name, database_name)
INTO {{PGDATABASE}}?load_ensembl_analysis_status
TARGET COLUMNS (task_name, database_name)
WITH fields terminated by ','

AFTER LOAD DO
$$
ALTER TABLE rnacen.load_ensembl_analysis_status SET (
    autovacuum_enabled = true,
    toast.autovacuum_enabled = true
);
$$,
$$
ANALYZE rnacen.load_ensembl_analysis_status;
$$
;
