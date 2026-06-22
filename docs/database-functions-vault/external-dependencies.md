# External dependencies

Functions that the release-load code calls from schemas **outside the ones this vault graphs**.

> **Note (2026-06-22):** as of the source-of-truth migration, all of these schemas
> (`rnc_logging`, `rnc_load_xref_incremental`, `database`, `upi`) **are now dumped into
> `database_functions/`** — their `.sql` files exist and deploy via the applier. They just
> don't have call-graph notes in this vault yet. The table below is still useful as the list of
> cross-schema call edges out of the release-load path.

| Schema | Function | Called from | Purpose (inferred) |
|---|---|---|---|
| `rnc_logging` | `log_release_start(dbid, release)` | [[rnc_update.load_release]] | Record start of a release load. |
| `rnc_logging` | `log_release_end(dbid, release, previous_release)` | [[rnc_update.load_release]] | Record completion / diff of a release load. |
| `rnc_load_xref_incremental` | `load_xref_incremental(previous_release, dbid)` | [[rnc_update.load_release]] | Xref load path for **incremental** (`release_type = 'I'`) releases; the full-release counterpart is [[rnc_load_xref.load_xref]]. |
| `database` | `set_current_release(dbid, release)` | [[rnc_update.mark_as_done]] | Point `rnc_database.current_release` at the finished release. |
| `upi` | `getUpi(seq_value)` | [[rnc_load_rna.load_md5_new_sequences_table]] | Format a `seq_upi` sequence value into a `URS……` UPI string. |

## Sequences (not functions)

Referenced via `nextval`/`currval`, so out of scope for a function call graph but worth noting:

- `seq_upi` — UPI allocation, used in [[rnc_load_rna.load_md5_new_sequences_table]].
- `xref_pk_seq` — xref primary keys, used in [[rnc_update.verify_xref_id_not_null]].
