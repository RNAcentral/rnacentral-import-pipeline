# rnc_load_rna.store_new_sequences

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`store_new_sequences__66805.sql`](../../../database_functions/rnc_load_rna/rnc_load_rna/store_new_sequences__66805.sql)

```sql
store_new_sequences() RETURNS void
```

Inserts the genuinely-new sequences into the main `rna` table, joining `load_retro_tmp` to
`load_md5_new_sequences` by md5 to pick up the allocated id/UPI. The final step of
[[rnc_load_rna.load_rna]].

- **Calls:** none
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `rna`, `load_retro_tmp`, `load_md5_new_sequences`
