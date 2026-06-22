# release.get_release_id

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_release_id__66758.sql`](../../../database_functions/release/release/get_release_id__66758.sql)

```sql
get_release_id(in_dbid bigint, in_release_date timestamp) RETURNS bigint
```

Returns the release id for a (database, release_date) pair. Swallows `no_data_found`.

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `rnc_release`
