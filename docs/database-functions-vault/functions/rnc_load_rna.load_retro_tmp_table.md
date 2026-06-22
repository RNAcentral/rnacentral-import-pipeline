# rnc_load_rna.load_retro_tmp_table

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`load_retro_tmp_table.sql`](../../../database_functions/rnc_load_rna/load_retro_tmp_table.sql)

```sql
load_retro_tmp_table(p_in_dbid bigint, p_in_load_release bigint) RETURNS void
```

Truncates and rebuilds `load_retro_tmp`. Selects distinct staged rows into a temp table,
left-joins them to existing `rna` by md5, and sets `comparable_prot_upi` to the existing
UPI when md5 + length + sequence all match (so re-uses of an existing sequence are detected).
The first step of [[rnc_load_rna.load_rna]].

- **Calls:** none
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `load_retro_tmp`, `load_rnacentral`, `rna`, `distinct_loaded_rows_tmp` (temp)
