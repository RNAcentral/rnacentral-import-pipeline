# release.get_load_release

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_load_release.sql`](../../../database_functions/release/get_load_release.sql)

```sql
get_load_release(in_dbid bigint) RETURNS bigint
```

Returns the release id for the database that is **not** done (`status <> 'D'`) — i.e. the
release currently being loaded. Swallows `no_data_found`.

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `rnc_release`
