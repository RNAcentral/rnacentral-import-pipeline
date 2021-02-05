COPY (
SELECT
    json_build_object(
        'id', todo.id,
        'urs_taxid', todo.urs_taxid,
        'urs_taxid', related.source_urs_taxid,
        'interacting_id', case when related.relationship_type = 'target_rna' then related.target_accession
                        else null
        'methods', related.methods
    )
FROM search_export_urs todo
JOIN rnc_related_sequences related
ON
  related.source_urs_taxid = todo.urs_taxid
WHERE
    related.relationship_type = 'target_rna'
ORDER BY todo.id
) TO STDOUT
