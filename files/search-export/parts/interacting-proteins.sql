COPY (
SELECT
    json_build_object(
        'id', related.source_urs_taxid,
        'interacting_protein_id', related.target_accession,
        'synonyms', protein.synonyms,
        'label', protein.label,
        'relationship', related.relationship_type,
        'methods', related.methods
    )
FROM rnc_related_sequences related
LEFT JOIN protein_info protein
ON
    protein.protein_accession = related.target_accession
WHERE
    related.relationship_type = 'target_protein'
) TO STDOUT
