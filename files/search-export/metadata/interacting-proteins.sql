COPY (
SELECT
    json_build_object(
        'rna_id', related.source_urs_taxid,
        'interacting_proteins', array_agg(
            json_build_object(
                'id', case
                        when related.relationship_type = 'target_protein'
                          then related.target_accession
                        else null
                      end,
                'synonyms', protein.synonyms,
                'label', protein.label,
                'relationship', related.relationship_type,
                'methods', related.methods
            )
        )
    )
FROM rnc_related_sequences related
LEFT JOIN protein_info protein
ON
    protein.protein_accession = related.target_accession
WHERE
    related.relationship_type = 'target_protein'
GROUP BY related.source_urs_taxid
) TO STDOUT
