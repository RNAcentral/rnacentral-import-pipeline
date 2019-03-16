LOAD CSV
FROM ALL FILENAMES MATCHING ~<matched_references.*csv$>
HAVING FIELDS (
    md5,
    pattern_group,
    pattern,
    matching_word,
    sentence,
    authors,
    location,
    title,
    pmid,
    doi
)
INTO {{PGDATABASE}}?load_matching_sentences
TARGET COLUMNS (
    md5,
    pattern_group,
    pattern,
    matching_word,
    sentence,
    authors,
    location,
    title,
    pmid,
    doi
)

WITH
    skip header = 0,
    fields escaped by double-quote,
    fields terminated by ','

BEFORE LOAD DO
$$
create table if not exists load_matching_sentences (
    pattern_group text,
    pattern text,
    matching_word text,
    sentence text,
    md5,
    authors,
    location,
    title,
    pmid,
    doi
);
$$

AFTER LOAD DO
$$
insert into rnc_references (
    md5,
    authors,
    location,
    title,
    pmid,
    doi
) (
    select distinct
        md5,
        authors,
        location,
        title,
        pmid,
        doi
    from load_matching_sentences
) ON CONFLICT (pmid) DO NOTHING;
$$,

$$
insert into rnc_matching_sentences (
    pattern_group,
    pattern,
    matching_word,
    sentence,
    reference_id
) (
    SELECT
        load.pattern_group,
        load.pattern,
        load.matching_word,
        load.sentence,
        refs.id
    from load_matching_sentences load
    join rnc_references refs on refs.pmid = load.pmid
);
$$,

$$
drop table load_matching_sentences
$$
;
