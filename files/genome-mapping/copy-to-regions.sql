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
select
	max(rna_id),
	region_id,
	max(chromosome),
	max(strand::int),
	min("start"),
	max(stop),
	assembly_id,
	true,
	max(mapping."identity"),
	'{}'
from rnc_genome_mapping mapping
group by mapping.region_id, mapping.assembly_id
) on conflict (region_name, assembly_id) do nothing
;
