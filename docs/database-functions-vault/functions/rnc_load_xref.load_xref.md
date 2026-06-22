# rnc_load_xref.load_xref

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`load_xref__66810.sql`](../../../database_functions/rnc_load_xref/rnc_load_xref/load_xref__66810.sql)

```sql
load_xref(p_previous_release bigint, p_in_dbid bigint) RETURNS void
```

Orchestrates a **full-release** xref load: build the pel tables, populate them with
kept/new/retired xrefs, then exchange them in as the live partitions and verify.

- **Calls (in order):**
  1. [[rnc_load_xref.prepare_pel_tables]]
  2. [[rnc_load_xref.populate_pel_tables1]]
  3. [[rnc_load_xref.populate_pel_tables2]]
  4. [[rnc_load_xref.load_upi_max_versions_table]]
  5. [[rnc_load_xref.load_max_versions_table]]
  6. [[rnc_load_xref.populate_pel_tables3]]
  7. [[rnc_load_xref.populate_pel_tables4]]
  8. [[rnc_load_xref.do_pel_exchange]]
- **Called by:** [[rnc_update.load_release]] (when `release_type = 'F'`)
- **Tables:** none directly

> **Perf change (2026-06-20):** the per-database `do_checks` call was removed.
> `do_checks` scans the *entire* `xref` table regardless of `dbid`, so calling it
> once per database meant re-checking all ~342M rows ~59 times per release. The
> pipeline (`release/run.py`) now patches this function at runtime and calls
> [[rnc_load_xref.do_checks]] **once after the whole per-database loop**.
