# rnc_load_xref.do_pel_exchange

Schema: [[rnc_load_xref]] · `SECURITY DEFINER`
Source: [`do_pel_exchange.sql`](../../../database_functions/rnc_load_xref/do_pel_exchange.sql)

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

> **Migrated to `DETACH/ATTACH PARTITION`** (declarative partitioning) from the old
> `INHERIT/NO INHERIT`. [[rnc_load_xref.revert_pel]] still uses the old scheme and is now stale
> (open item #1).
>
> **Perf change (2026-06-22, [[follow-ups]] §3b):** the `upi → rna` FKs (`fk4` on each
> partition) are added `NOT VALID` here — the data satisfies them by construction, so the
> in-transaction full-partition scan + probes into the 35 GB `rna(upi)` index was pure overhead
> (~hours for dbid 1). `release/run.py` runs `VALIDATE CONSTRAINT` after the per-DB loop
> (`ShareUpdateExclusiveLock`, off the critical path). Still rebuilds 14 indexes per exchange.
