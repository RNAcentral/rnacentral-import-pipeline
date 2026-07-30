# rnc_load_xref.populate_pel_tables1

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`populate_pel_tables1.sql`](../../../database_functions/rnc_load_xref/populate_pel_tables1.sql)

```sql
populate_pel_tables1(v_previous_release bigint) RETURNS void
```

For existing non-deleted xrefs that match a staged row, decides whether each stays active
or is retired (CTE `l_data`), then routes rows to `xref_pel_not_deleted` (deleted `'N'`) or
`xref_pel_deleted` (deleted `'Y'`). Handles the "unchanged / re-seen" entries.

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `xref`, `load_retro_tmp`, `xref_pel_deleted`, `xref_pel_not_deleted`

> **Perf change (2026-06-20):** the `xref` join (`x.dbid = l.in_dbid`) is a join
> key, so the planner can't prune the LIST(dbid) partitions and Appends across all
> ~59. `release/run.py` patches this function at runtime to read the single dbid
> from `load_retro_tmp` into a local `v_dbid` and add `and x.dbid = v_dbid` — a
> redundant predicate (load_retro_tmp is single-dbid) that enables partition
> pruning (verified `Subplans Removed: 58`). Same fix applied to
> [[rnc_load_xref.populate_pel_tables2]] and [[rnc_load_xref.load_upi_max_versions_table]].
