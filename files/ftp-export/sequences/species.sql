COPY (
SELECT DISTINCT
    pre.taxid,
    tax.name
FROM rnc_rna_precomputed pre
JOIN rnc_taxonomy tax ON tax.id = pre.taxid
WHERE
    pre.is_active = true
    AND pre.taxid IS NOT NULL
ORDER BY tax.name
) TO STDOUT CSV
