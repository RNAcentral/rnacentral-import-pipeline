COPY (
SELECT
    hits.urs,
    hits.rfam_model_id,
    score,
    e_value,
    sequence_start,
    sequence_stop,
    model_start,
    model_stop,
    models.long_name
FROM rfam_model_hits hits
JOIN rnc_rna_precomputed pre on pre.urs = hits.urs and taxid is null
JOIN rfam_models models ON models.rfam_model_id = hits.rfam_model_id
WHERE
    pre.is_active = true
    AND hits.sequence_stop > hits.sequence_start
ORDER BY hits.urs, hits.sequence_start, hits.rfam_model_id
) TO STDOUT
