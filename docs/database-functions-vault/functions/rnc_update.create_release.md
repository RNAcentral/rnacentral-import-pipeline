# rnc_update.create_release

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`create_release__66853.sql`](../../../database_functions/rnc_update/rnc_update/create_release__66853.sql)

```sql
create_release(p_in_dbid bigint, p_release_type character) RETURNS void
```

Inserts a single new `rnc_release` row for a database: id = `count(*)+1`, today's date,
status `'L'`, userstamp `auto`.

- **Calls:** none
- **Called by:** [[rnc_update.prepare_releases]]
- **Tables:** `rnc_release`
