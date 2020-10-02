COPY (
SELECT
    json_build_object(
        'id', related.source_urs_taxid,
        'interacting_rna_id', related.target_accession,
        'urs', related.target_urs_taxid,
        'methods', related.methods
    )
FROM rnc_related_sequences related
WHERE
    related.relationship_type = 'target_rna'
) TO STDOUT
