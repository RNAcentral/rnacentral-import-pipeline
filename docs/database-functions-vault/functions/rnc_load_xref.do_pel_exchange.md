# rnc_load_xref.do_pel_exchange

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`do_pel_exchange__66807.sql`](../../../database_functions/rnc_load_xref/rnc_load_xref/do_pel_exchange__66807.sql)

```sql
do_pel_exchange(p_in_db_id bigint) RETURNS void
```

The partition swap. Renames the current `xref_p<dbid>_(not_)deleted` partitions (and their
indexes/FKs) to `_old`, detaches them, renames the freshly built `xref_pel_*` tables into
place, adds check constraints / indexes / FKs, reattaches them to the right partition parent,
and analyzes. [[rnc_load_xref.revert_pel]] is its inverse.

- **Calls:** none
- **Called by:** [[rnc_load_xref.load_xref]]
- **Tables:** `xref` partitions (`xref_p<dbid>_deleted`, `xref_p<dbid>_not_deleted`), `xref_pel_*`; catalogs `pg_constraint`, `pg_inherits`
