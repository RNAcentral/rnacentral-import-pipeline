LOAD CSV
FROM ALL FILENAMES MATCHING ~<matched_references.*csv$>
HAVING FIELDS (
    md5,
    accession,
    authors,
    location,
    title,
    pmid,
    doi,
    word,
    sentence,
    source
)
INTO {{PGDATABASE}}?rnc_matching_sentences
TARGET COLUMNS (
    md5,
    accession,
    authors,
    location,
    title,
    pmid,
    doi,
    word,
    sentence,
    source
)

WITH
    truncate,
    drop indexes,
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','
;
