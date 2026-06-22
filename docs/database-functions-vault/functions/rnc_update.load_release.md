# rnc_update.load_release

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`load_release__66854.sql`](../../../database_functions/rnc_update/rnc_update/load_release__66854.sql)

```sql
load_release(p_in_dbid bigint, p_in_load_release bigint) RETURNS void
```

The heart of a release load. Logs start, loads sequences, then loads xrefs down the **full**
or **incremental** path depending on the release type, marks the release done, and logs end.

- **Calls (in order):**
  - `rnc_logging.log_release_start` — see [[external-dependencies]]
  - [[rnc_load_rna.load_rna]]
  - [[release.get_release_type]]
  - [[release.get_previous_release]]
  - [[rnc_load_xref.load_xref]] *(when type = `F`)*
  - `rnc_load_xref_incremental.load_xref_incremental` *(when type = `I`)* — see [[external-dependencies]]
  - [[rnc_update.mark_as_done]]
  - `rnc_logging.log_release_end` — see [[external-dependencies]]
- **Called by:** [[rnc_update.new_update_release]]
- **Tables:** none directly
