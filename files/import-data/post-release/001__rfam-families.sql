\timing

BEGIN;

INSERT INTO rfam_models (
    rfam_model_id,
    short_name,
    long_name,
    description,
    seed_count,
    full_count,
    length,
    is_suppressed,
    rfam_clan_id,
    domain,
    rna_type,
    rfam_rna_type,
    so_rna_type
) (
SELECT
    rfam_model_id,
    short_name,
    long_name,
    description,
    seed_count,
    full_count,
    length,
    is_suppressed,
    rfam_clan_id,
    domain,
    rna_type,
    rfam_rna_type,
    so_rna_type
FROM load_rfam_models
)
ON CONFLICT (rfam_model_id) DO UPDATE SET
    short_name = EXCLUDED.short_name,
    long_name = EXCLUDED.long_name,
    description = EXCLUDED.description,
    seed_count = EXCLUDED.seed_count,
    full_count = EXCLUDED.full_count,
    length = EXCLUDED.length,
    is_suppressed = EXCLUDED.is_suppressed,
    rfam_clan_id = EXCLUDED.rfam_clan_id,
    domain = EXCLUDED.domain,
    rna_type = EXCLUDED.rna_type,
    rfam_rna_type = EXCLUDED.rfam_rna_type
;

DROP TABLE load_rfam_models;

COMMIT;
