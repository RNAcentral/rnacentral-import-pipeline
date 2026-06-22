# rnc_update.move_staging_data

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`move_staging_data__66856.sql`](../../../database_functions/rnc_update/rnc_update/move_staging_data__66856.sql)

```sql
move_staging_data(p_in_dbid bigint) RETURNS void
```

Truncates `load_rnacentral` and repopulates it with the distinct rows of
`load_rnacentral_all` belonging to this database (matched via `rnc_database.descr`). Narrows
the multi-database staging table down to the one being loaded.

- **Calls:** none
- **Called by:** [[rnc_update.new_update_release]]
- **Tables:** `load_rnacentral`, `load_rnacentral_all`, `rnc_database`
