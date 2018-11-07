INSERT INTO ref_pubmed (
    ref_pubmed_id,
    authors,
    location,
    title,
    doi
) (
SELECT DISTINCT
    ref_pubmed_id,
    authors,
    location,
    title,
    doi
FROM load_ref_pubmed
)
ON CONFLICT (ref_pubmed_id) DO UPDATE
SET
    authors = EXCLUDED.authors,
    location = EXCLUDED.location,
    title = EXCLUDED.title,
    doi = EXCLUDED.title
;

DROP TABLE load_ref_pubmed;
