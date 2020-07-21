LOAD CSV
FROM ALL FILENAMES MATCHING ~<traveler-data.*.csv$>
HAVING FIELDS (
    urs,
    model,
    secondary_structure,
    layout,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    stk
) INTO {{PGDATABASE}}?rnc_secondary_structure_layout
TARGET COLUMNS (
    urs,
    model,
    secondary_structure,
    layout,
    overlap_count,
    basepair_count,
    model_start,
    model_stop,
    sequence_start,
    sequence_stop,
    sequence_coverage,
    stk
)

WITH
    batch rows = 300,
    batch concurrency = 3,
    FIELDS ESCAPED BY double-quote,
    FIELDS TERMINATED BY ','
;
