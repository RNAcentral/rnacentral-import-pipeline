# rnc_load_xref.populate_pel_tables4

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`populate_pel_tables4.sql`](../../../database_functions/rnc_load_xref/populate_pel_tables4.sql)

```sql
populate_pel_tables4(p_in_dbid bigint, v_previous_release bigint) RETURNS void
```

Writes the **deleted/retired** xrefs into `xref_pel_deleted`: existing rows already marked
deleted (`UNION ALL`) plus active rows whose accession is absent from this load. Sets `last`
to the previous release for entries being retired now.

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `xref_pel_deleted`, `xref`, `load_retro_tmp`
