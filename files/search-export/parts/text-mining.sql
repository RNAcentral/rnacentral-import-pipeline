\o /dev/null
CREATE TEMPORARY TABLE search_export_publication_counts
AS
SELECT
  UPPER(d.primary_id) as urs,
  SUM(j.hit_count) as publication_count
FROM litscan_job j
JOIN litscan_database d
ON d.job_id = j.job_id
WHERE j.hit_count > 0 AND d.name='rnacentral'
GROUP BY d.primary_id;

\o

COPY (
  SELECT
    json_build_object(
      'id', todo.id,
      'urs_taxid', todo.urs_taxid,
      'publication_count', counts.publication_count
    )
    FROM search_export_urs todo
    JOIN search_export_publication_counts counts
    ON
      todo.urs_taxid = counts.urs
    ORDER by todo.id
) TO STDOUT
