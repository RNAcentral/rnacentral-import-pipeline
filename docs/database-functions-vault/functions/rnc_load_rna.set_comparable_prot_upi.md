# rnc_load_rna.set_comparable_prot_upi

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`set_comparable_prot_upi.sql`](../../../database_functions/rnc_load_rna/set_comparable_prot_upi.sql)

```sql
set_comparable_prot_upi() RETURNS void
```

Fills `comparable_prot_upi` in `load_retro_tmp` for the new sequences, copying the freshly
allocated UPI from `load_md5_new_sequences` where it was still NULL. Then rebuilds the
`load_retro_tmp$ac_dbid_upi` index. After this every staged row has a UPI.

- **Calls:** none
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `load_retro_tmp`, `load_md5_new_sequences`

> **Perf change (2026-06-22):** rewritten from a per-row correlated subquery to a single
> `UPDATE load_retro_tmp l SET comparable_prot_upi = n.prot_upi FROM load_md5_new_sequences n
> WHERE n.in_md5 = l.in_md5 AND l.comparable_prot_upi IS NULL` (hash join; equivalent because
> `in_md5` is unique in `load_md5_new_sequences`, and non-matching rows stay NULL either way).
> The supporting `load_md5_new_sequences$in_md5` index (added 2026-06-20 when this was still a
> correlated subquery) is kept — the planner may hash-join instead, but it's a harmless safety net.
