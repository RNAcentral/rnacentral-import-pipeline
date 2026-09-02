LOAD CSV
FROM ALL FILENAMES MATCHING ~<feedback.*tsv$>
HAVING FIELDS (
    upi_taxid,
    status,
    result,
    assembly_id
)
INTO {{PGDATABASE}}?load_overlaps
TARGET COLUMNS (
    upi_taxid,
    status,
    result,
    assembly_id
)

WITH
    fields terminated by ' '

AFTER LOAD DO
$$
insert into rnc_feedback_overlap (
    urs_taxid,
    overlaps_with,
    overlapping_urs_taxids,
    no_overlaps_with,
    assembly_id
) (
select
    load.upi_taxid,
    coalesce(
        array_agg(distinct load.result) FILTER (WHERE load.status = 'overlap'),
        '{}'::text[]
    ),
    coalesce(
        array_agg(distinct load.result) FILTER (WHERE load.status = 'overlapping_id'),
        '{}'::text[]
    ),
    coalesce(
        array_agg(distinct load.result) FILTER (WHERE load.status = 'no_overlap'),
        '{}'::text[]
    ),
    load.assembly_id
from load_overlaps load
join rnc_rna_precomputed pre on pre.urs_taxid = load.upi_taxid
group by load.upi_taxid
) ON CONFLICT (urs_taxid, assembly_id) DO UPDATE
SET
    overlaps_with = excluded.overlaps_with,
    overlapping_urs_taxids = excluded.overlapping_urs_taxids,
    no_overlaps_with = excluded.no_overlaps_with
;
$$,
$$
TRUNCATE load_overlaps;
$$
;
