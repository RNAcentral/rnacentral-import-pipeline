INSERT INTO rfam_clans (
    rfam_clan_id,
    name,
    description,
    family_count
) (
SELECT
    rfam_clan_id,
    name,
    description,
    family_count
FROM load_rfam_clans
)
ON CONFLICT (rfam_clan_id) DO UPDATE SET
    name = EXCLUDED.name,
    description = EXCLUDED.description,
    family_count = EXCLUDED.family_count
;

DROP TABLE load_rfam_clans;
