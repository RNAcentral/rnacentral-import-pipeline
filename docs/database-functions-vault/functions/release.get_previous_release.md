# release.get_previous_release

Schema: [[release]] · `STABLE SECURITY DEFINER`
Source: [`get_previous_release__66757.sql`](../../../database_functions/release/release/get_previous_release__66757.sql)

```sql
get_previous_release(in_dbid bigint, in_release_id bigint) RETURNS bigint
```

Returns the largest release id less than `in_release_id` for the database.

- **Calls:** none
- **Called by:** [[release.get_retired_count]], [[rnc_update.load_release]]
- **Tables:** `rnc_release`
