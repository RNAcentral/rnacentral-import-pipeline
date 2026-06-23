COPY (
SELECT DISTINCT(ontology_term_id)
FROM ontology_terms
WHERE ontology = 'SO'
) TO STDOUT
