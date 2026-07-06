# release.get_retired_count

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_retired_count.sql`](../../../database_functions/release/get_retired_count.sql)

```sql
get_retired_count(in_dbid bigint, in_release bigint) RETURNS integer
```

Counts xrefs retired during `in_release`: those whose `last` equals the previous release id.

- **Calls:** [[release.get_previous_release]]
- **Called by:** application code (not from within this folder)
- **Tables:** `xref`
