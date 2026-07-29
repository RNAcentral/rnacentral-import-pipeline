# release.get_release_status

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_release_status.sql`](../../../database_functions/release/get_release_status.sql)

```sql
get_release_status(in_release_id bigint) RETURNS character
```

Returns the `status` of a release. **Raises** `no_data_found` if the release does not exist.

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `rnc_release`
