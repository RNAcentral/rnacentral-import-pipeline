# rnc_load_rna.load_md5_stats_table

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`load_md5_stats_table.sql`](../../../database_functions/rnc_load_rna/load_md5_stats_table.sql)

```sql
load_md5_stats_table() RETURNS void
```

Rebuilds `load_md5_stats`: for each md5 among the *new* staged rows (those in
`load_retro_tmp` with `comparable_prot_upi IS NULL`), counts rows and distinct short/long
sequences. Feeds collision detection and new-sequence allocation.

- **Calls:** none
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `load_md5_stats`, `load_retro_tmp`
