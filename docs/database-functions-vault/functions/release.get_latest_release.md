# release.get_latest_release

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_latest_release.sql`](../../../database_functions/release/get_latest_release.sql)

```sql
get_latest_release(in_dbid bigint) RETURNS bigint
```

Returns the maximum release id for the database. Swallows `no_data_found`.

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `rnc_release`
