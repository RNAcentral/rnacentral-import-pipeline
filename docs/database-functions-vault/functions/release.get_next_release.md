# release.get_next_release

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_next_release.sql`](../../../database_functions/release/get_next_release.sql)

```sql
get_next_release(in_dbid bigint, in_release_id bigint) RETURNS bigint
```

Returns the smallest release id greater than `in_release_id` for the database.

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `rnc_release`
