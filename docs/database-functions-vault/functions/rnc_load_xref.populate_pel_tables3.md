# rnc_load_xref.populate_pel_tables3

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`populate_pel_tables3__66813.sql`](../../../database_functions/rnc_load_xref/rnc_load_xref/populate_pel_tables3__66813.sql)

```sql
populate_pel_tables3(p_in_dbid bigint) RETURNS void
```

Inserts brand-new and re-appearing accessions into `xref_pel_not_deleted`. Computes
`version_i` from `load_max_versions` / `load_upi_max_versions`: `1` for new, unchanged when
the UPI is the same, `+1` when the UPI changed.

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `xref_pel_not_deleted`, `load_upi_max_versions`, `load_retro_tmp`, `load_max_versions`
