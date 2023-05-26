LOAD CSV
FROM ALL FILENAMES MATCHING ~<merged_metadata>
HAVING FIELDS (
    job_id,
    name,
    primary_id
) INTO {{PGDB_EMBASSY_USER}}?litscan_database
TARGET COLUMNS (
    job_id,
    name,
    primary_id
)

WITH
    batch rows = 1000000,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY '|'

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS litscan_database;
$$,
$$
CREATE TABLE litscan_database (
    id SERIAL PRIMARY KEY,
    job_id VARCHAR,
    name VARCHAR,
    primary_id VARCHAR,
    FOREIGN KEY (job_id) REFERENCES litscan_job(job_id) ON UPDATE CASCADE ON DELETE CASCADE,
    FOREIGN KEY (primary_id) REFERENCES litscan_job(job_id) ON UPDATE CASCADE ON DELETE CASCADE,
    CONSTRAINT name_job UNIQUE (job_id, name, primary_id)
);
$$

AFTER LOAD DO
$$
CREATE INDEX ON litscan_database (name);
$$,
$$
GRANT SELECT ON litscan_database TO rnacen;
$$
;
