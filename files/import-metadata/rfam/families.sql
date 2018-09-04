select
    family.rfam_acc id,
    family.rfam_id as name,
    family.description as pretty_name,
    group_concat(concat('SO:', dbs.db_link)) as so_terms,
    family.type as rna_type,
    family.comment as description,
    family.num_seed as seed_count,
    family.num_full as full_count,
    group_concat(mem.clan_acc) as clan_id,
    family.clen as length
from family
left join database_link dbs on dbs.rfam_acc = family.rfam_acc
left join clan_membership mem on mem.rfam_acc = family.rfam_acc
where
    dbs.db_id = 'SO'
group by family.rfam_acc

