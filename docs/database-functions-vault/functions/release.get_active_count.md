# release.get_active_count

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_active_count__66751.sql`](../../../database_functions/release/release/get_active_count__66751.sql)

```sql
get_active_count(in_dbid bigint, in_release bigint) RETURNS integer
```

Counts xrefs for `in_dbid` that are active at `in_release` (created on/before it and either
not deleted or last-seen at/after it).

- **Calls:** none
- **Called by:** application code (not from within this folder)
- **Tables:** `xref`
