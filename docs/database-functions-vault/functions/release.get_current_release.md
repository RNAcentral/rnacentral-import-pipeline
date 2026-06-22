# release.get_current_release

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_current_release__66753.sql`](../../../database_functions/release/release/get_current_release__66753.sql)

```sql
get_current_release(in_dbid bigint) RETURNS bigint
```

Returns `rnc_database.current_release` for the database. Swallows `no_data_found`.

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `rnc_database`
