COPY (
select
  database_name, task_name
from ensembl_import_tracking tracking
where
  tracking.was_imported = true
) TO STDOUT CSV
