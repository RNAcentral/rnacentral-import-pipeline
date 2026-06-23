# rnc_load_rna.load_rna

Schema: [[rnc_load_rna]] · `STABLE SECURITY DEFINER`
Source: [`load_rna.sql`](../../../database_functions/rnc_load_rna/load_rna.sql)

```sql
load_rna(p_in_dbid bigint, p_in_load_release bigint) RETURNS void
```

Orchestrates loading new sequences into `rna`, running the six steps below in order.

- **Calls (in order):**
  1. [[rnc_load_rna.load_retro_tmp_table]]
  2. [[rnc_load_rna.load_md5_stats_table]]
  3. [[rnc_load_rna.load_md5_collisions_table]]
  4. [[rnc_load_rna.load_md5_new_sequences_table]]
  5. [[rnc_load_rna.set_comparable_prot_upi]]
  6. [[rnc_load_rna.store_new_sequences]]
- **Called by:** [[rnc_update.load_release]]
- **Tables:** none directly
