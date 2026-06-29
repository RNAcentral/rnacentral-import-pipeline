# rnc_update.verify_xref_id_not_null

Schema: [[rnc_update]] · `SECURITY DEFINER`
Source: [`verify_xref_id_not_null.sql`](../../../database_functions/rnc_update/verify_xref_id_not_null.sql)

```sql
verify_xref_id_not_null() RETURNS void
```

Backfills any `xref` rows with a null `id`, assigning the next value of `xref_pk_seq`.

- **Calls:** none
- **Called by:** nothing active — only referenced in a **commented-out** line of
  [[rnc_load_xref.do_checks]]
- **Tables:** `xref` · **Sequence:** `xref_pk_seq`
