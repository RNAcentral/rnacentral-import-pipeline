COPY (
SELECT
    json_build_object(
        'rna_id', related.source_urs_taxid,
        'interacting_rnas', array_agg(
            json_build_object(
                'id', case
                        when related.relationship_type = 'target_rna'
                          then related.target_accession
                        else null
                      end,
                'urs', related.target_urs_taxid,
                'methods', related.methods
            )
        )
    )
FROM rnc_related_sequences related
WHERE
    related.relationship_type = 'target_rna'
GROUP BY related.source_urs_taxid
) TO STDOUT
