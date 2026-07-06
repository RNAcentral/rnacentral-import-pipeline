# rnc_update.update_rnc_accessions

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`update_rnc_accessions.sql`](../../../database_functions/rnc_update/update_rnc_accessions.sql)

```sql
update_rnc_accessions() RETURNS void
```

Standalone loader. Upserts `rnc_accessions` from `load_rnc_accessions`
(`INSERT … ON CONFLICT (accession) DO UPDATE` across the descriptive columns), then analyzes
the table. Sets `application_name` and logs EXPLAIN output.

> Note: selects `t2.so_term` into the `rna_type` column and `t2.rna_type` is used in the
> conflict-update — worth keeping in mind if column mappings ever change.

- **Calls:** none
- **Called by:** pipeline / operator (entry point)
- **Tables:** `rnc_accessions`, `load_rnc_accessions`
