# rnc_load_xref.do_checks

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`do_checks__66806.sql`](../../../database_functions/rnc_load_xref/rnc_load_xref/do_checks__66806.sql)

```sql
do_checks(p_in_db_id bigint) RETURNS void
```

Post-exchange sanity check. Builds/refreshes the materialized view `xref_pk_not_unique`
(ids appearing more than once across the inheriting tables) and raises `unique_violation`
if any duplicates exist. Final step of [[rnc_load_xref.load_xref]].

> Contains a **commented-out** call to [[rnc_update.verify_xref_id_not_null]] — kept here as
> a documented (inactive) edge.

- **Calls:** none active (commented: `rnc_update.verify_xref_id_not_null`)
- **Called by:** the pipeline (`release/run.py`) **once after the per-database loop**.
  No longer called from [[rnc_load_xref.load_xref]] (removed 2026-06-20 — see that note).
- **Tables / objects:** `xref`, materialized view `xref_pk_not_unique`

> The `p_in_db_id` argument is **ignored** — the query (`select id, count(*) from
> xref group by id having count(*) > 1`) is global. `run.py` calls it as
> `do_checks(NULL::bigint)`. This is a candidate for removal: if the run after
> this change shows `xref_pk_not_unique` is always empty, the check (likely a
> leftover debugging aid from the inheritance-partitioning era) can be dropped.
