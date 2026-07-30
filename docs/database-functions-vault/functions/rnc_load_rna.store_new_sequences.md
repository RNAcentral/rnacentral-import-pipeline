# rnc_load_rna.store_new_sequences

Schema: [[rnc_load_rna]] · `SECURITY DEFINER`
Source: [`store_new_sequences.sql`](../../../database_functions/rnc_load_rna/store_new_sequences.sql)

```sql
store_new_sequences() RETURNS void
```

Inserts the genuinely-new sequences into the main `rna` table, joining `load_retro_tmp` to
`load_md5_new_sequences` by md5 to pick up the allocated id/UPI. The final step of
[[rnc_load_rna.load_rna]].

- **Calls:** none
- **Called by:** [[rnc_load_rna.load_rna]]
- **Tables:** `rna`, `load_retro_tmp`, `load_md5_new_sequences`

> **Perf change (2026-06-22):** this was the `load_rna` killer (see [[follow-ups]] §6). The
> original `SELECT DISTINCT … seq_long` deduped on the full long-sequence text, so the sort
> detoasted/compared multi-KB sequences per row — a multi-day CPU spin that emitted zero rows
> (the INSERT can't start until the sort finishes). Changed to `DISTINCT ON (in_md5)`, which
> sorts only the 32-char md5. Equivalent: an md5 only reaches `load_md5_new_sequences` when
> `load_md5_stats` saw a single distinct sequence for it, so all rows sharing an md5 carry the
> same sequence (and would otherwise collide on the `rna.id` PK).
