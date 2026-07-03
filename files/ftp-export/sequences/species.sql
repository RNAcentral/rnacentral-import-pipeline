COPY (
SELECT
    pre.taxid,
    replace(tax.name, ',', '') AS tax_name,
    replace(tax.lineage, ',', ';') AS lineage
FROM rnc_rna_precomputed pre
JOIN rnc_taxonomy tax ON tax.id = pre.taxid
WHERE
    pre.is_active = true
    AND pre.taxid IS NOT NULL
GROUP BY pre.taxid, tax.name, tax.lineage
HAVING count(*) > 100
ORDER BY tax_name
) TO STDOUT CSV
