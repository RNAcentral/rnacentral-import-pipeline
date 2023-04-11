COPY(
  SELECT
    article.pmcid,
    organism.organism
  FROM litscan_article article
  JOIN litscan_load_organism organism
  ON
    article.pmid = organism.pmid
  WHERE article.pmid IS NOT NULL
  ORDER BY article.pmcid
) TO STDOUT (FORMAT CSV)
