COPY (
SELECT DISTINCT
	CASE rfam.rfam_model_id IS NULL
	WHEN true THEN 'rrna'
	WHEN false THEN 'rfam'
	END AS model_type,
	ont."name" AS rna_type
FROM rnc_secondary_structure_layout_models models 
LEFT JOIN rfam_models rfam ON rfam.rfam_model_id = models.model_name
JOIN ontology_terms ont ON ont.ontology_term_id = models.so_term_id
) TO STDOUT CSV
