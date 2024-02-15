BEGIN;

-- Create a table for the loaded assemblies
CREATE TEMP TABLE mapped_assemblies AS SELECT DISTINCT assembly_id FROM load_genome_mapping;
CREATE INDEX ix_mapped_assemblies ON mapped_assemblies (assembly_id);

-- Create a temp table for all possible assemblies to track changes in counts
-- from
CREATE TEMP TABLE all_assemblies AS
SELECT DISTINCT assembly_id FROM rnc_genome_mapping;
CREATE UNIQUE INDEX un_all_assemblies__assembly_id ON all_assemblies (assembly_id);

INSERT INTO all_assemblies (assembly_id) (
  SELECT DISTINCT assembly_id FROM load_genome_mapping
) ON CONFLICT (assembly_id) DO NOTHING;

-- Index load_mapping for faster access
CREATE INDEX ix_load_genome_mapping__assembly_id ON load_genome_mapping (assembly_id);

DECLARE current_version int;
SET current_version = select coalesce(max(update_version), 0) from rnc_genome_mapping_stats;

INSERT INTO rnc_genome_mapping_stats (
  assembly_id,
  update_version,
  previous_count,
  current_count,
  new_regions,
  deleted_regions,
  common_regions
) (
  SELECT
    assem.assembly_id,
    :current_version,
    (
      select
        count(*)
      from rnc_sequence_regions regions
      where
        regions.assembly_id = assem.assembly_id
        and regions.was_mapped = true
    ) as previous_count,
    (
      select
        count(distinct mapping.region_name)
      from load_genome_mapping mapping
      where
        mapping.assembly_id = assem.assembly_id
    ) as current_count,
    (
      select
        count(distinct mapping.region_name)
      from load_genome_mapping mapping
      left join rnc_sequence_regions regions on regions.region_name = mapping.region_name
      where
        mapping.assembly_id = assem.assembly_id
        and regions.id is null
    ) as new_regions,
    (
      select
        count(distinct regions.id)
      from rnc_sequence_regions regions
      left join load_genome_mapping mapping
      on
        regions.region_name = mapping.region_name
        and mapping.assembly_id = assem.assembly_id
      where
        regions.assembly_id = assem.assembly_id
        and regions.was_mapped = true
        and mapping.region_name is null
    ) as deleted_regions,
    (
      select
        count(distinct regions.region_name)
      from load_genome_mapping mapping
      join rnc_sequence_regions regions on regions.region_name = mapping.region_name
      where
        mapping.assembly_id = assem.assembly_id
        and regions.was_mapped = true
    ) as common_regions
  from all_assemblies assem
  group by assem.assembly_id
)
;

DELETE FROM rnc_sequence_regions regions
USING mapped_assemblies load
WHERE
    regions.assembly_id = load.assembly_id
    AND regions.was_mapped = true
;

INSERT INTO rnc_sequence_regions (
    urs_taxid,
    region_name,
    chromosome,
    strand,
    region_start,
    region_stop,
    assembly_id,
    exon_count,
    was_mapped,
    identity,
    providing_databases
) (
SELECT
    max(load.urs_taxid),
    load.region_name,
    max(load.chromosome),
    max(load.strand),
    min(load.exon_start),
    max(load.exon_stop),
    load.assembly_id,
    max(load.exon_count),
    true,
    max(load.identity),
    '{}'::text[]
FROM load_genome_mapping load
JOIN ensembl_assembly ensembl on ensembl.assembly_id = load.assembly_id
GROUP BY load.region_name, load.assembly_id
) ON CONFLICT (MD5(region_name), assembly_id) DO NOTHING
;

INSERT INTO rnc_sequence_exons (
    region_id,
    exon_start,
    exon_stop
) (
SELECT
    regions.id,
    load.exon_start,
    load.exon_stop
FROM load_genome_mapping load
JOIN rnc_sequence_regions regions ON regions.region_name = load.region_name
JOIN ensembl_assembly ensembl ON ensembl.assembly_id = load.assembly_id
) ON CONFLICT (region_id, exon_start, exon_stop) DO NOTHING
;

DROP TABLE load_genome_mapping;

COMMIT;
