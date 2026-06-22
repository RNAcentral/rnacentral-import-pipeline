# rnc_update.mark_as_done

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`mark_as_done__66855.sql`](../../../database_functions/rnc_update/rnc_update/mark_as_done__66855.sql)

```sql
mark_as_done(p_in_dbid bigint, p_in_load_release bigint) RETURNS void
```

Finalises a release: sets its status to `'D'` and points the database's current release at it.

- **Calls:**
  - [[release.set_release_status]]
  - `database.set_current_release` — see [[external-dependencies]]
- **Called by:** [[rnc_update.load_release]]
- **Tables:** none directly
