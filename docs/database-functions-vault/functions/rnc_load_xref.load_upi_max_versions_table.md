# rnc_load_xref.load_upi_max_versions_table

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`load_upi_max_versions_table.sql`](../../../database_functions/rnc_load_xref/load_upi_max_versions_table.sql)

```sql
load_upi_max_versions_table(p_in_dbid bigint) RETURNS void
```

Rebuilds `load_upi_max_versions`: for staged accessions that have **no** active xref yet,
finds the max previous `version_i` per (ac, dbid) along with the UPI it belonged to. Used to
continue version numbering for re-appearing accessions.

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `load_upi_max_versions`, `xref`, `load_retro_tmp`

> **Perf change (2026-06-20):** patched at runtime (via `release/run.py`) so the
> `xref` predicates reference the `p_in_dbid` argument directly
> (`PREVIOUS_XREF.DBID = p_in_dbid`, and `X.DBID = p_in_dbid` in the `NOT EXISTS`)
> instead of the `l.in_dbid` join column. The `WHERE` already enforces
> `L.IN_DBID = p_in_dbid`, so this is exact — and it lets the planner prune the
> LIST(dbid) partitions. See [[rnc_load_xref.populate_pel_tables1]].
