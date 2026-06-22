# rnc_load_rna.set_comparable_prot_upi

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`set_comparable_prot_upi__66804.sql`](../../../database_functions/rnc_load_rna/rnc_load_rna/set_comparable_prot_upi__66804.sql)

```sql
set_comparable_prot_upi() RETURNS void
```

Fills `comparable_prot_upi` in `load_retro_tmp` for the new sequences, copying the freshly
allocated UPI from `load_md5_new_sequences` where it was still NULL. Then rebuilds the
`load_retro_tmp$ac_dbid_upi` index. After this every staged row has a UPI.

- **Calls:** none
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `load_retro_tmp`, `load_md5_new_sequences`

> **Perf change (2026-06-20):** the correlated subquery probes
> `load_md5_new_sequences` by `in_md5` once per `load_retro_tmp` row, and that
> table had no index (per-row seq scan). `release/run.py` now creates
> `load_md5_new_sequences$in_md5` once at setup (the table is only `TRUNCATE`d, not
> dropped, per load, so the index persists across the per-database loop).
