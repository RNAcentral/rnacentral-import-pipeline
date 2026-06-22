# rnc_load_rna.load_md5_collisions_table

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`load_md5_collisions_table.sql`](../../../database_functions/rnc_load_rna/load_md5_collisions_table.sql)

```sql
load_md5_collisions_table() RETURNS void
```

Rebuilds `load_md5_collisions` from `load_md5_stats`, keeping md5s that map to more than one
distinct short or long sequence (i.e. genuine hash collisions worth flagging).

- **Calls:** none
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `load_md5_collisions`, `load_md5_stats`
