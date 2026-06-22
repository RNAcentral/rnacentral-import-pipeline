# release.get_release_type

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_release_type__66760.sql`](../../../database_functions/release/release/get_release_type__66760.sql)

```sql
get_release_type(in_release_id bigint) RETURNS character
```

Returns the `release_type` of a release (e.g. `F` full / `I` incremental). **Raises**
`no_data_found` if missing. The returned type decides which xref path
[[rnc_update.load_release]] takes.

- **Calls:** none
- **Called by:** [[rnc_update.load_release]]
- **Tables:** `rnc_release`
