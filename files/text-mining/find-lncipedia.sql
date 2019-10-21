COPY (
SELECT
    acc.gene
FROM xref
JOIN rnc_accessions acc 
ON 
    acc.accession = xref.ac
WHERE
    deleted = 'N'
    AND dbid = 20
) TO STDOUT
