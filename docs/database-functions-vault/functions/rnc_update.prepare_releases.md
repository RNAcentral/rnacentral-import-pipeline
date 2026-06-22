# rnc_update.prepare_releases

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`prepare_releases.sql`](../../../database_functions/rnc_update/prepare_releases.sql)

```sql
prepare_releases(p_release_type character) RETURNS void
```

Seeds the release table before a load. If any release is already in `'L'` (loading) state it
returns early; otherwise it creates one new release per database present in the staging table.

- **Calls:** [[rnc_update.create_release]] (once per staged database)
- **Called by:** pipeline / operator (entry point)
- **Tables:** `rnc_release`, `load_rnacentral_all`, `rnc_database`
