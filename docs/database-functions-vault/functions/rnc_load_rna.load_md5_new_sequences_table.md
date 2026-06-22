# rnc_load_rna.load_md5_new_sequences_table

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`load_md5_new_sequences_table.sql`](../../../database_functions/rnc_load_rna/load_md5_new_sequences_table.sql)

```sql
load_md5_new_sequences_table() RETURNS void
```

Rebuilds `load_md5_new_sequences`: for each genuinely-new md5 (no collision in
`load_md5_stats`), allocates a new id from the `seq_upi` sequence and formats it into a
`URS……` UPI. This is where new sequences get their permanent identifier.

- **Calls:** `upi.getUpi(...)` — see [[external-dependencies]]
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `load_md5_new_sequences`, `load_md5_stats` · **Sequence:** `seq_upi`
