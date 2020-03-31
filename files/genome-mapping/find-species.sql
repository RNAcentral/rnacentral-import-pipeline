COPY (
SELECT
  species.ensembl_url,
  species.assembly_id,
  species.taxid,
  species.division
FROM :species_to_map species
where
  exists(select 1 from :tablename to_map where to_map.assembly_id = species.assembly_id)
) TO STDOUT CSV;
