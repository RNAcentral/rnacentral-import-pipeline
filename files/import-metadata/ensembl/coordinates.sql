select
	seq_region.name,
	coord.name as coordinate_system,
	meta.meta_value,
	attrib_type.code as attrib_name,
	seq_region_attrib.value as attrib_value
from seq_region
join coord_system coord on coord.coord_system_id = seq_region.coord_system_id
join meta on meta.species_id = coord.species_id
join seq_region_attrib on seq_region_attrib.seq_region_id = seq_region.seq_region_id
join attrib_type on attrib_type.attrib_type_id = seq_region_attrib.attrib_type_id
where
	coord.name = 'chromosome'
	and meta_key = 'assembly.default'
order by seq_region.name
;
