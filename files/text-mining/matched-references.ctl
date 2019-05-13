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
INTO {{PGDATABASE}}?load_rnc_text_mining
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
CREATE TABLE IF NOT EXISTS load_rnc_text_mining (
    pattern_group text,
    pattern text,
    matching_word text,
    sentence text,
    md5 text,
    authors text,
    location text,
    title text,
    pmid text,
    doi text
);
$$

AFTER LOAD DO
$$
INSERT INTO rnc_references (
    md5,
    authors,
    location,
    title,
    pmid,
    doi
) (
    SELECT distinct
        md5,
        authors,
        location,
        title,
        pmid,
        doi
    FROM load_rnc_text_mining
) ON CONFLICT (md5) DO NOTHING;
$$,

$$
INSERT INTO rnc_text_mining_patterns (
  pattern_group,
  pattern_name
) (
SELECT
  pattern_group, 
  pattern
FROM load_rnc_text_mining
) ON CONFLICT (pattern_group, pattern_name) DO NOTHING;
$$,

$$
INSERT INTO rnc_text_mining_sentences (
  sentence,
  reference_id
) (
SELECT
  load.sentence, 
  refs.id
FROM load_rnc_text_mining load
JOIN rnc_references refs
ON
  refs.pmid = load.pmid
) ON CONFLICT (md5(sentence), reference_id) DO NOTHING;
$$,

$$
INSERT INTO rnc_text_mining_matches (
    matching_word,
    pattern_id,
    sentence_id
) (
    SELECT
        load.matching_word,
        patt.id,
        sent.id
    FROM load_rnc_text_mining load
    JOIN rnc_text_mining_sentences sent
    ON
      sent.sentence = load.sentence
    JOIN rnc_references refs 
    ON 
      refs.id = sent.reference_id
      AND refs.pmid = load.pmid
    JOIN rnc_text_mining_patterns patt
    ON
      patt.pattern_group = load.pattern_group
      AND patt.pattern_name = load.pattern
) ON CONFLICT (matching_word, pattern_id, sentence_id) DO NOTHING;
$$
;
