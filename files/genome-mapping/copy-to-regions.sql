insert into rnc_sequence_regions (
	urs_taxid,
	region_name,
	chromosome,
	strand,
	region_start,
	region_stop,
	assembly_id,
	was_mapped,
	identity,
	providing_databases
) (
select distinct
	max(rna_id),
	region_id,
	max(chromosome),
	max(strand::int),
	min("start"),
	max(stop),
	assembly_id,
	true,
	max(mapping."identity"),
	'{}'::text[]
from rnc_genome_mapping mapping
group by mapping.region_id, mapping.assembly_id
) on conflict (region_name, assembly_id) do nothing
;

insert into rnc_sequence_exons (
  region_id,
  exon_start,
  exon_stop
) (
select distinct
  regions.id,
  mapping.start,
  mapping.stop
from rnc_genome_mapping mapping
join rnc_sequence_regions regions on mapping.region_id = regions.region_name and regions.assembly_id = mapping.assembly_id
)
;
