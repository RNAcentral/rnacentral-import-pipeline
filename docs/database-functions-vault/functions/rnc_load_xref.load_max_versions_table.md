# rnc_load_xref.load_max_versions_table

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`load_max_versions_table.sql`](../../../database_functions/rnc_load_xref/load_max_versions_table.sql)

```sql
load_max_versions_table() RETURNS void
```

Collapses `load_upi_max_versions` to one row per (ac, dbid) holding the maximum `version_i`,
storing it in `load_max_versions`. Consumed by [[rnc_load_xref.populate_pel_tables3]].

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `load_max_versions`, `load_upi_max_versions`
