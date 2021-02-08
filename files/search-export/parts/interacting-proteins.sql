COPY (
SELECT
    json_build_object(
        'id', todo.id,
        'urs_taxid', todo.urs_taxid,
        'interacting_protein_id', case
                        when related.relationship_type = 'target_protein'
                          then related.target_accession
                        else null
                      end,
        'synonyms', protein.synonyms,
        'label', protein.label,
        'relationship', related.relationship_type,
        'methods', related.methods
    )
FROM search_export_urs todo
join rnc_related_sequences related
ON
  related.source_urs_taxid = todo.urs_taxid
LEFT JOIN protein_info protein
ON
    protein.protein_accession = related.target_accession
WHERE
    related.relationship_type = 'target_protein'
ORDER BY todo.id
) TO STDOUT
