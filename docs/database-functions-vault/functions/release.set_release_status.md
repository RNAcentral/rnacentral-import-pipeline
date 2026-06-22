# release.set_release_status

Schema: [[release]] · `SECURITY DEFINER`
Source: [`set_release_status.sql`](../../../database_functions/release/set_release_status.sql)

```sql
set_release_status(in_release_id bigint, release_status character) RETURNS void
```

Updates a release's `status` column. Used to flip a release to `'D'` (done) at the end of a load.

- **Calls:** none
- **Called by:** [[rnc_update.mark_as_done]]
- **Tables:** `rnc_release`
