COPY (
SELECT
    pre.taxid,
    replace(tax.name, ',', '') AS tax_name
FROM rnc_rna_precomputed pre
JOIN rnc_taxonomy tax ON tax.id = pre.taxid
WHERE
    pre.is_active = true
    AND pre.taxid IS NOT NULL
GROUP BY pre.taxid, tax.name
HAVING count(*) > 10
ORDER BY tax_name
) TO STDOUT CSV
