# rnc_load_xref.populate_pel_tables2

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`populate_pel_tables2.sql`](../../../database_functions/rnc_load_xref/populate_pel_tables2.sql)

```sql
populate_pel_tables2() RETURNS void
```

Inserts the **updated** active xrefs into `xref_pel_not_deleted`: matched staged rows that
were *not* simply re-seen (the negated condition that distinguishes this from
[[rnc_load_xref.populate_pel_tables1]]). Bumps `version_i` when the UPI changed.

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `xref_pel_not_deleted`, `load_retro_tmp`, `xref`

> **Perf change (2026-06-20):** patched at runtime (via `release/run.py`) to add a
> `and x.dbid = v_dbid` partition-pruning predicate on the `xref` join. See
> [[rnc_load_xref.populate_pel_tables1]] for the full rationale.
