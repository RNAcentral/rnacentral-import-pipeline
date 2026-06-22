# rnc_load_xref.revert_pel

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`revert_pel__66816.sql`](../../../database_functions/rnc_load_xref/rnc_load_xref/revert_pel__66816.sql)

```sql
revert_pel(p_in_db_id bigint) RETURNS void
```

Manual rollback / inverse of [[rnc_load_xref.do_pel_exchange]]. Swaps the `inherit` clauses
back, drops the just-installed partitions, and renames the `_old` partitions, indexes and FK
constraints back to current. **Not** part of the [[rnc_load_xref.load_xref]] flow — run by
hand to undo a bad exchange.

- **Calls:** none
- **Called by:** nothing in this folder (manual)
- **Tables:** `xref` partitions (`xref_p<dbid>_(not_)deleted` and their `_old` variants)
