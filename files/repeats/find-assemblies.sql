COPY (
SELECT
  species.ensembl_url,
  species.assembly_id,
  species.division
FROM ensembl_assembly species
WHERE
  exists(select 1 from rnc_sequence_regions reg where reg.assembly_id = species.assembly_id)
  and species.division != 'EnsemblFungi'
) TO STDOUT CSV;

