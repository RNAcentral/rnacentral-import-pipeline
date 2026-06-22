# rnc_update.update_literature_references

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`update_literature_references__66859.sql`](../../../database_functions/rnc_update/rnc_update/update_literature_references__66859.sql)

```sql
update_literature_references() RETURNS void
```

Standalone loader. Inserts new rows into `rnc_references` (distinct refs from
`load_rnc_references` not already present by md5), inserts `(accession, reference_id)` links
into `rnc_reference_map` with `ON CONFLICT DO NOTHING`, then truncates the staging table.
Sets `application_name` and logs EXPLAIN output along the way.

- **Calls:** none
- **Called by:** pipeline / operator (entry point)
- **Tables:** `rnc_references`, `rnc_reference_map`, `load_rnc_references`
