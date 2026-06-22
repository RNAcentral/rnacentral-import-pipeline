# release.get_checkpoint

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_checkpoint.sql`](../../../database_functions/release/get_checkpoint.sql)

```sql
get_checkpoint() RETURNS bigint
```

Returns the highest release id whose `status = 'D'` (the last fully-done release).

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `rnc_release`
