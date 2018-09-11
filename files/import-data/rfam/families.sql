SELECT
    family.rfam_acc id,
    family.rfam_id AS name,
    family.description AS pretty_name,
    group_concat(concat('SO:', dbs.db_link)) AS so_terms,
    family.type AS rna_type,
    family.comment AS description,
    family.num_seed AS seed_count,
    family.num_full AS full_count,
    group_concat(mem.clan_acc) AS clan_id,
    family.clen AS length
FROM family
LEFT JOIN clan_membership mem ON mem.rfam_acc = family.rfam_acc
WHERE
    dbs.db_id = 'SO'
GROUP BY family.rfam_acc
