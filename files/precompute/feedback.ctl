LOAD CSV
FROM 'merged.tsv'
HAVING FIELDS (
    upi_taxid,
    status,
    result
)
INTO {{PGDATABASE}}?load_overlaps
TARGET COLUMNS (
    upi_taxid,
    status,
    result
)

WITH
    fields terminated by ' '

BEFORE LOAD DO
$$
drop table if exists load_overlaps;
$$,
$$
create table load_overlaps (
    upi_taxid text,
    status text,
    result text
);
$$

AFTER LOAD DO
$$
insert into rnc_feedback_overlap (
    upi_taxid,
    overlaps_with,
    overlapping_upis,
    no_overlaps_with
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
    )
from load_overlaps load
join rnc_rna_precomputed pre on pre.id = load.upi_taxid
group by load.upi_taxid
) ON CONFLICT (upi_taxid) DO UPDATE
SET
    overlaps_with = excluded.overlaps_with,
    overlapping_upis = excluded.overlapping_upis,
    no_overlaps_with = excluded.no_overlaps_with
;
$$
;
