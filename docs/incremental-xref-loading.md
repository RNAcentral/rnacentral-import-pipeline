# Xref loading: full rebuild vs delta — design & correctness notes

**Branch:** `improvement/load-only-new-data`
**Status:** two release types: `'F'` (full rebuild) and `'D'` (delta, for a
manifest-tracked database — currently ENA only). Selection is per-database
(default: auto) via `release.get_load_release_type`: bootstrap, and any database
without an import manifest, get `'F'`; a manifest-tracked database gets `'D'`.
The remaining work is **staging-scale validation + benchmarking** (step 4) before
delta is trusted on production volumes for a given database. `rnac release run
--force-full` reverts the release to full for every database at any time;
`--force_full_import` on the pipeline does that *and* forces a full parse, which
is what a delta-parsed database needs (see docs/incremental-parsing.md).

## Problem

Loading is slow because every release rewrites *every* xref for each database,
even rows that are unchanged. We want to load only rows that are **new** or have
**genuinely changed** (including going inactive), leaving unchanged rows in place.

## How loading works today

Two layers:

1. **Staging load** — [`bin/split-and-load`](../bin/split-and-load) runs pgloader
   with a `.ctl` file (e.g. [`accessions.ctl`](../files/import-data/load/accessions.ctl))
   that **truncates and bulk-loads the full CSV** into `load_*` tables every run.
   Cost is set by what the parsers emit — a full dump per database, except ENA's
   delta-parsed export (see docs/incremental-parsing.md).

2. **Release merge** — `rnac release run`
   ([`run.py`](../rnacentral_pipeline/rnacentral/release/run.py)) moves staging →
   production. Two sub-parts:
   - **Sequences (`rna`)** — already incremental.
     [`store_new_sequences`](../database_functions/rnc_load_rna/store_new_sequences.sql)
     inserts only sequences whose md5 is new. Unchanged sequences are never rewritten.
   - **Xrefs (`xref`)** — either a full rebuild (`'F'`) or a delta (`'D'`), chosen
     per database.

### Full rebuild (`'F'`)

A full release runs Partition Exchange Loading (PEL):

- `populate_pel_tables1..4` build brand-new partition tables
  (`xref_pel_not_deleted` / `xref_pel_deleted`) containing **all** xrefs for the
  database — active *and* the entire history of deleted rows.
- [`do_pel_exchange`](../database_functions/rnc_load_xref/do_pel_exchange.sql)
  swaps them in, replacing the whole partition.

So a full release rewrites, for each database:
- every unchanged active row (only `last` actually changes), **and**
- **every historical `deleted='Y'` row, carried forward untouched**
  ([`populate_pel_tables4`](../database_functions/rnc_load_xref/populate_pel_tables4.sql),
  first `UNION ALL` arm).

That deleted-history carry-forward grows without bound and is pure waste every
release — unavoidable for `'F'`, but exactly what `'D'`'s row-by-row updates skip
for a database whose input really is a small delta.

## Delta (`'D'`): correctness rules

`'D'` mutates `xref` in place instead of rebuilding it, so its steps must
reproduce the same end state a full rebuild would produce for the rows they
touch. Ground truth = the state the full/PEL path produces. Per `(ac, dbid)`,
given the new load row `L` (with `comparable_prot_upi`, `version`, `taxid`,
`in_load_release`) and the previous release `P`:

| Case | Condition | Required end state |
|---|---|---|
| 1. Unchanged active | active row `X`, `X.urs = L.urs`, `X.version = L.version` | keep row; `last = in_load_release`; `deleted='N'`; `taxid = COALESCE(L.taxid, X.taxid)`; `version_i`, `created` unchanged |
| 2. Changed active | active `X`, `X.ac=L.ac` but urs/version differ | retire `X` in place (`deleted='Y'`, `last = P`); **insert** new active row: `urs=L.urs`, `version_i = X.version_i (+1 if urs changed)`, `created=last=in_load_release`, `deleted='N'` |
| 3. New sequence for existing accession | accession has active rows but none with `urs=L.urs` (PEL "Gap A") | insert active row, `version_i = max(version_i)+1` |
| 4. Brand-new accession / previously fully deleted | no active row for `(ac,dbid)` | insert active row, `version_i = 1` if none ever, else `max+1` (keep if same urs) |
| 5. Explicitly deleted | accession named in the parser's `load_deletions` list | retire in place: `deleted='Y'`, `last = COALESCE(P, X.last)` |
| 6. Absent from the load, not explicitly deleted | active `X` whose `(ac, urs)` is simply not in this delta | **leave active** — absence means "unchanged", not "gone" |
| 7. Already deleted | `X.deleted='Y'` | leave untouched (do **not** rewrite — this is the win) |

Field rules:
- **`version_i`** — monotonic per `(ac, dbid)`. New → 1; same `comparable_prot_upi`
  → unchanged; changed → `max+1` (matches `populate_pel_tables2`).
- **`deleted`** — `'N'` active / `'Y'` retired; at most one active row per generation.
- **`taxid`** — refresh = `COALESCE(incoming, existing)`; otherwise unchanged.
- **`last`** — surviving active row bumped to `in_load_release` (marks "seen this
  release"); retired row set to previous release.
- **`created`** — set once at insert, never changed.

## Choosing full vs delta

- **First load of a database** (no prior completed release / no partition) →
  **FULL** (bootstraps the PEL partition tables).
- **A database with an import manifest** (its parser emits only new/changed
  records — see docs/incremental-parsing.md) → **DELTA**, except **HGNC**, whose
  auto-DELTA is deliberately paused for now (see docs/incremental-parsing.md) —
  it runs FULL despite having a real manifest.
- **Every other subsequent load** → **FULL**. A row-by-row pass over a full-size
  input costs no less than a rebuild and does more work per row (see risk 2
  below), so there's no in-between mode for a database whose parser doesn't
  actually produce a delta.
- **Escape hatch** — `rnac release run --force-full` forces `'F'` for every
  database at any time (schema change, suspected drift, or after changing the
  merge logic).

Implementation: `release.get_load_release_type(dbid)` returns `'F'` when there is
no prior release or no import manifest for the dbid, `'D'` otherwise — with an
HGNC-specific override ahead of that check, pinning it to `'F'` regardless.
`rnc_update.prepare_releases` has an `'A'` (auto) mode that applies that choice
per database; an explicit `'F'`/`'D'` still forces a type for every database.

## Risks / open questions

1. **Parity.** `'D'`'s row-by-row steps must yield the identical active set a full
   rebuild would, for the rows they touch. Mitigation:
   `tests/xref-incremental-parity/delta.sh` exercises the case table above end to
   end (unchanged-and-absent stays active, explicit deletion retires, version
   change replaces, new accession inserts) and diffs against the expected state.
2. **When most rows change**, in-place UPDATE/INSERT (index churn, dead tuples,
   autovacuum) can be *slower* than build-fresh-and-swap — `'D'` only wins for a
   genuinely small delta, which is why it's gated behind an import manifest
   rather than applied by default.
3. **Reader impact.** PEL swap is invisible until exchange; `'D'` mutates the live
   partition (row locks, bloat) while the site reads it. Consider vacuum/timing.
4. **Transactionality.** PEL swaps atomically. `'D'`'s mutations must stay in one
   transaction so a mid-release failure can't half-update a partition.
   `new_update_release` is a single function call (one txn) — verify under the
   `run.py` autocommit-per-statement driver.
5. **Downstream.** `populate_precompute` and `do_checks` read the active set / PK
   uniqueness — both must hold identically after a delta load.
6. **fk4 (`urs → rna`).** Full validates fk4 per partition post-swap; `'D'`
   inserts into an already-valid partition (checked at insert). Safe as long as
   `store_new_sequences` runs first (it does).

## Findings that make the in-place approach safe

- **`xref` is declaratively partitioned** (`dbid` LIST → `deleted` LIST; see
  [`do_pel_exchange`](../database_functions/rnc_load_xref/do_pel_exchange.sql):
  `attach partition … for values in ('Y'/'N')`). Therefore `INSERT INTO xref`
  auto-routes to the right partition, and an `UPDATE` that flips `deleted='N'→'Y'`
  **moves the row** from `_not_deleted` to `_deleted` (PG 11+ row movement). The
  in-place approach needs no manual partition juggling.
- **Full releases already reassign every `id`.** The PEL inserts
  ([`populate_pel_tables*`](../database_functions/rnc_load_xref/)) omit `id`, so a
  rebuilt partition gets fresh ids for every carried-forward row. `xref.id` is thus
  **not stable across full releases**; `'D'` keeps ids stable, which is strictly
  better. ⇒ **Parity must be compared on the business key**
  `(ac, dbid, version, version_i, urs, deleted, last, taxid)`, never on `id`. The
  only uniqueness invariant [`do_checks`](../database_functions/rnc_load_xref/do_checks.sql)
  enforces is `id` uniqueness, trivially preserved (new rows take the default id).
- **Derivation strategy:** each row-by-row step is the **in-place analogue of its
  `populate_pel_tables*` counterpart** (same predicates and value expressions,
  retargeted from the `xref_pel_*` build tables to `rnacen.xref`). This makes
  parity hold *by construction*. Because `create_release` sets the new release id
  to `count(*)+1` (strictly greater than every prior id), **all pre-existing rows
  satisfy `last < in_load_release`**, so a `last < in_load_release` guard on the
  UPDATE steps cleanly excludes rows this same release just inserted — with no
  false exclusions.

### Step-by-step mapping (`'D'` ⇐ PEL)

Run order (helper snapshot first, then inserts, then in-place updates):

1. `load_upi_max_versions_table(dbid)` + `load_max_versions_table()` — snapshot of
   the *original* xref (must precede all mutations).
2. **A1** `incremental_new_versions` ⇐ `populate_pel_tables2`: INSERT a new active
   generation for accessions whose sequence version changed (same urs).
3. **A2/A3** `incremental_new_accessions` ⇐ `populate_pel_tables3` (main + Gap A):
   INSERT active rows for brand-new / previously-fully-deleted accessions and for a
   new sequence variant of an already-active accession.
4. **B1** `incremental_refresh` ⇐ `populate_pel_tables1` (deleted='N' branch):
   in-place UPDATE of matched-unchanged rows (`last`, `taxid`).
5. **B2** `incremental_retire_changed` ⇐ `populate_pel_tables1` (deleted='Y'
   branch): in-place UPDATE retiring matched-but-changed rows.
6. **B3** `incremental_retire_explicit` ⇐ the parser's `load_deletions` list:
   in-place UPDATE retiring only the accessions explicitly named as gone. Absence
   alone never retires under `'D'` — that's the entire performance win, and the
   reason `'D'` needs a manifest-tracked parser to be correct at all.

`populate_pel_tables4`'s *first* arm (carry forward all historical `deleted='Y'`
rows) has **deliberately no `'D'` analogue** — leaving those rows untouched is
the rest of the performance win.

## Proposed sequence of work

1. ✅ **Row-by-row steps built** — `incremental_new_versions`,
   `incremental_new_accessions`, `incremental_refresh`, `incremental_retire_changed`
   (each the in-place analogue of a `populate_pel_tables*` step) and
   `incremental_retire_explicit` (retires only the parser's explicit deletion
   list), orchestrated by `load_xref_delta`.
2. ✅ **Delta parity validated** — `tests/xref-incremental-parity/delta.sh` seeds a
   known xref state, runs `load_xref_delta`, and asserts every case in the table
   above: absence stays active, explicit deletion retires, version change
   replaces, new accession inserts, already-deleted history is untouched.
3. ✅ **Selection logic** — `release.get_load_release_type(dbid)` returns `'F'` for
   a database with no prior release (bootstrap) or no import manifest, `'D'` for
   one with a manifest. `rnc_update.prepare_releases` gains an `'A'` (auto) mode
   that applies that choice per database; an explicit `'F'`/`'D'` still forces a
   type. `run.run()` calls `prepare_releases('A')` by default, and
   `rnac release run --force-full` forces `'F'` for every database. Tested by:
   - `tests/rnacentral/release/run_test.py` — `run()` issues `prepare_releases('A')`
     by default and `prepare_releases('F')` under `--force-full`.
   - `tests/xref-incremental-parity/selection.sh` — DB-level assertions that
     `get_load_release_type` and `prepare_releases('A'|'F')` assign the right types.
4. ⬜ **Staging validation + benchmark** — run delta against a real staging copy
   of ENA; measure wall-clock vs full; watch bloat/autovacuum on the mutated live
   partition. **Do this before relying on `'D'` for further databases**; until
   then, `--force-full` keeps the old full-rebuild behaviour available at any time.
5. ⬜ **Roll out delta parsing to more databases** — extend the manifest/delta
   parsing approach (docs/incremental-parsing.md) beyond ENA once staging
   validation on ENA holds up.

### Validating the harness locally

    PGHOST=127.0.0.1 PGPORT=5433 PGUSER=postgres ./tests/xref-incremental-parity/delta.sh
    PGHOST=127.0.0.1 PGPORT=5433 PGUSER=postgres ./tests/xref-incremental-parity/selection.sh

(any Postgres you can create databases in.)

## Out of scope (noted)

The staging pgloader step still truncate-reloads the full CSV for every database
except ENA, because their parsers emit a full dump. Making that incremental would
require each parser to diff against DB state — larger change, separate effort.
