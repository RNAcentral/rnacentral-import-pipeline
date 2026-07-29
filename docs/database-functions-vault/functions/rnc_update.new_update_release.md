# rnc_update.new_update_release

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`new_update_release.sql`](../../../database_functions/rnc_update/new_update_release.sql)

```sql
new_update_release(p_in_dbid bigint, p_in_release_id bigint,
                   p_release_type character DEFAULT 'f') RETURNS void
```

**Main entry point** for running a release update for one database. Moves staging data, then
loads the named release.

- **Calls:**
  - [[rnc_update.move_staging_data]]
  - [[rnc_update.load_release]]
- **Called by:** pipeline / operator (top-level)
- **Tables:** declares a cursor over `rnc_release` (status `'L'`)
