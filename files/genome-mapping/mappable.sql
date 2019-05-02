COPY (
SELECT
  ensembl_url,
  assembly_id,
  taxid,
  division
FROM ensembl_assembly
WHERE
  division != 'EnsemblProtists'
  AND division != 'EnsemblFungi'
  AND EXISTS(SELECT 1 
    FROM rnc_rna_precomputed pre 
    WHERE 
      pre.taxid = ensembl_assembly.taxid 
      AND pre.has_coordinates IS false 
      AND pre.is_active IS true
  )
) TO STDOUT CSV;
