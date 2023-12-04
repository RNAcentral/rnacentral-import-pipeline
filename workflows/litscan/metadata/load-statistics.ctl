LOAD CSV
FROM ALL FILENAMES MATCHING ~<statistics.csv>
HAVING FIELDS (
    searched_ids,
    articles,
    ids_in_use,
    urs,
    expert_db
) INTO {{PGDATABASE}}?litscan_statistics
TARGET COLUMNS (
    searched_ids,
    articles,
    ids_in_use,
    urs,
    expert_db
)

WITH
    FIELDS TERMINATED BY ','

BEFORE LOAD DO
$$
DROP TABLE IF EXISTS litscan_statistics;
$$,
$$
CREATE TABLE litscan_statistics (
    id SERIAL PRIMARY KEY,
    searched_ids INT,
    articles INT,
    ids_in_use INT,
    urs INT,
    expert_db INT
);
$$

AFTER LOAD DO
$$
GRANT SELECT ON litscan_statistics TO rnacen;
$$
;
