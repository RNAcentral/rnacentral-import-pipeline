INSERT INTO ensembl_import_tracking (
    database_name,
    task_name,
    was_imported
) (
SELECT
    database_name,
    task_name,
    true
FROM load_ensembl_analysis_status
) ON CONFLICT (database_name, task_name) DO UPDATE
SET
    was_imported = excluded.was_imported
;

DROP TABLE load_ensembl_analysis_status;
