# Incremental xref loading — design & correctness notes

**Branch:** `improvement/load-only-new-data`
**Status:** incremental xref functions written and **parity-validated**, and the
**per-database FULL/INCREMENTAL selection is now wired in** (default: auto). The
remaining work is **staging-scale validation + benchmarking** (step 5) before this
is trusted on production volumes. `rnac release run --force-full` reverts to the
old all-full behaviour at any time.

## Problem

Loading is slow because every release rewrites *every* xref for each database,
even rows that are unchanged. We want to load only rows that are **new** or have
**genuinely changed** (including going inactive), leaving unchanged rows in place.

## How loading works today

Two layers:

1. **Staging load** — [`bin/split-and-load`](../bin/split-and-load) runs pgloader
   with a `.ctl` file (e.g. [`accessions.ctl`](../files/import-data/load/accessions.ctl))
   that **truncates and bulk-loads the full CSV** into `load_*` tables every run.
   Cost is set by what the parsers emit (a full dump per database). Comparatively
   cheap; **not** the main bottleneck and **not** changed by this proposal.

2. **Release merge** — `rnac release run`
   ([`run.py`](../rnacentral_pipeline/rnacentral/release/run.py)) moves staging →
   production. Two sub-parts:
   - **Sequences (`rna`)** — **already incremental**.
     [`store_new_sequences`](../database_functions/rnc_load_rna/store_new_sequences.sql)
     inserts only sequences whose md5 is new. Unchanged sequences are never rewritten. ✅
   - **Xrefs (`xref`)** — **fully rebuilt every release**. This is the bottleneck.

### Why xrefs are the bottleneck

[`run.py:107`](../rnacentral_pipeline/rnacentral/release/run.py#L107) hard-codes
`prepare_releases('F')` — **every** release is type `'F'` (full). A full release
runs Partition Exchange Loading (PEL):

- `populate_pel_tables1..4` build brand-new partition tables
  (`xref_pel_not_deleted` / `xref_pel_deleted`) containing **all** xrefs for the
  database — active *and* the entire history of deleted rows.
- [`do_pel_exchange`](../database_functions/rnc_load_xref/do_pel_exchange.sql)
  swaps them in, replacing the whole partition.

So every release rewrites, for each database:
- every unchanged active row (only `last` actually changes), **and**
- **every historical `deleted='Y'` row, carried forward untouched**
  ([`populate_pel_tables4`](../database_functions/rnc_load_xref/populate_pel_tables4.sql),
  first `UNION ALL` arm).

That deleted-history carry-forward grows without bound and is pure waste every
release. It is exactly what an incremental load avoids.

## The existing (disabled, incomplete) incremental path

Release type `'I'` already exists:
[`load_xref_incremental`](../database_functions/rnc_load_xref_incremental/load_xref_incremental.sql)
→ `incremental1/2/3`. It is **never invoked** (always `'F'`) and has two problems:

1. **It is Oracle SQL and will not run on Postgres.** e.g. in
   [`incremental1`](../database_functions/rnc_load_xref_incremental/incremental1.sql):
   `UPDATE XREF U SET (U.LAST, U.DELETED, u.taxid) = (SELECT …)` — Postgres forbids
   table-qualified LHS and uses different multi-column-UPDATE syntax; plus `NVL`
   (→ `COALESCE`), `/*+ PARALLEL(...) */` hints (drop), and `X."VERSION" = NULL`
   (→ `IS NULL`).

2. **It is functionally incomplete — it never deactivates disappeared xrefs.**
   `incremental1` only touches rows that are *matched* by the new load.
   The full path handles "previously active, not in the new load → `deleted='Y'`"
   in `populate_pel_tables4`; the incremental path has **no equivalent**. An
   accession dropped from a source would stay `active` forever. This is precisely
   the "become inactive" behaviour we want, so a new `incremental4` step is
   required — it does not exist yet.

## Correctness rules the incremental path MUST reproduce

Ground truth = the state the full/PEL path produces. Per `(ac, dbid)`, given the
new load row `L` (with `comparable_prot_upi`, `version`, `taxid`, `in_load_release`)
and the previous release `P`:

| Case | Condition | Required end state |
|---|---|---|
| 1. Unchanged active | active row `X`, `X.upi = L.upi`, `X.version = L.version` | keep row; `last = in_load_release`; `deleted='N'`; `taxid = COALESCE(L.taxid, X.taxid)`; `version_i`, `created` unchanged |
| 2. Changed active | active `X`, `X.ac=L.ac` but upi/version differ | retire `X` in place (`deleted='Y'`, `last = P`); **insert** new active row: `upi=L.upi`, `version_i = X.version_i (+1 if upi changed)`, `created=last=in_load_release`, `deleted='N'` |
| 3. New sequence for existing accession | accession has active rows but none with `upi=L.upi` (PEL "Gap A") | insert active row, `version_i = max(version_i)+1` |
| 4. Brand-new accession / previously fully deleted | no active row for `(ac,dbid)` | insert active row, `version_i = 1` if none ever, else `max+1` (keep if same upi) |
| 5. Disappeared | active `X` whose `(ac[,upi])` is **not** in new load | retire in place: `deleted='Y'`, `last = COALESCE(P, X.last)` — **missing from current incremental** |
| 6. Already deleted | `X.deleted='Y'` | leave untouched (do **not** rewrite — this is the win) |

Field rules:
- **`version_i`** — monotonic per `(ac, dbid)`. New → 1; same `comparable_prot_upi`
  → unchanged; changed → `max+1`.
  ⚠️ **Divergence to reconcile:** `incremental2` bumps `version_i` on *version*
  change too, whereas `populate_pel_tables2` bumps only on *upi* change. Pick one
  (match PEL) and apply consistently.
- **`deleted`** — `'N'` active / `'Y'` retired; at most one active row per generation.
- **`taxid`** — refresh = `COALESCE(incoming, existing)`; otherwise unchanged.
- **`last`** — surviving active row bumped to `in_load_release` (marks "seen this
  release"); retired row set to previous release.
- **`created`** — set once at insert, never changed.

## Choosing full vs incremental

- **First load of a database** (no prior completed release / no partition) → must be
  **FULL** (bootstraps the PEL partition tables).
- **Subsequent loads** → **INCREMENTAL**.
- **Escape hatch** → `rnc_release.force_load` already exists; honour it to force a
  FULL rebuild (schema change, suspected drift, or after changing the merge logic).

Implementation: replace the hard-coded `prepare_releases('F')` with per-database
type selection — `'I'` when a prior completed release exists for the dbid and
`force_load='N'`, else `'F'`. `create_release`/`get_release_type`/
`get_previous_release` already support this; only the *selection* is hard-coded.

## Risks / open questions

1. **Parity.** Incremental must yield the identical active set to full. Two known
   divergences (case 5 missing, `version_i` bump rule). Mitigation: a diff harness
   that runs full vs incremental on the same input and compares the `active` xref
   set — must be built before switching the default.
2. **Oracle → Postgres rewrite** of `incremental1/2/3` + new `incremental4`.
3. **When most rows change**, in-place UPDATE/INSERT (index churn, dead tuples,
   autovacuum) can be *slower* than build-fresh-and-swap. Incremental wins for the
   common small-delta case; consider a per-database threshold or measurement.
4. **Reader impact.** PEL swap is invisible until exchange; incremental mutates the
   live partition (row locks, bloat) while the site reads it. Consider vacuum/timing.
5. **Transactionality.** PEL swaps atomically. Incremental must wrap each database's
   mutations in one transaction so a mid-release failure can't half-update a
   partition. `new_update_release` is a single function call (one txn) — verify
   under the `run.py` autocommit-per-statement driver.
6. **Downstream.** `populate_precompute` and `do_checks` read the active set / PK
   uniqueness — both must hold identically after incremental.
7. **fk4 (`upi → rna`).** Full validates fk4 per partition post-swap; incremental
   inserts into an already-valid partition (checked at insert). Safe as long as
   `store_new_sequences` runs first (it does).

## Findings that make the in-place approach safe (verified while reviewing)

- **`xref` is declaratively partitioned** (`dbid` LIST → `deleted` LIST; see
  [`do_pel_exchange`](../database_functions/rnc_load_xref/do_pel_exchange.sql):
  `attach partition … for values in ('Y'/'N')`). Therefore `INSERT INTO xref`
  auto-routes to the right partition, and an `UPDATE` that flips `deleted='N'→'Y'`
  **moves the row** from `_not_deleted` to `_deleted` (PG 11+ row movement). The
  in-place incremental approach needs no manual partition juggling.
- **Full releases already reassign every `id`.** The PEL inserts
  ([`populate_pel_tables*`](../database_functions/rnc_load_xref/)) omit `id`, so a
  rebuilt partition gets fresh ids for every carried-forward row. `xref.id` is thus
  **not stable across full releases today**; incremental keeps ids stable, which is
  strictly better. ⇒ **Parity must be compared on the business key**
  `(ac, dbid, version, version_i, upi, deleted, last, taxid)`, never on `id`. The
  only uniqueness invariant [`do_checks`](../database_functions/rnc_load_xref/do_checks.sql)
  enforces is `id` uniqueness, trivially preserved (new rows take the default id).
- **`check_function_bodies` is effectively off** in the deploy path
  ([`functions.apply`](../rnacentral_pipeline/rnacentral/release/functions.py) applies
  every `database_functions/*/*.sql` in one transaction). That is why the broken
  Oracle `incremental*.sql` deploy without error and only fail at *runtime* — which
  never happens because no `'I'` release is ever created. Replacing them is safe,
  and callee-before-caller apply order is not a concern.
- **Derivation strategy chosen:** rather than porting the Oracle originals, each
  incremental step is written as the **in-place analogue of its
  `populate_pel_tables*` counterpart** (same predicates and value expressions,
  retargeted from the `xref_pel_*` build tables to `rnacen.xref`). This makes
  parity hold *by construction*. Because `create_release` sets the new release id
  to `count(*)+1` (strictly greater than every prior id), **all pre-existing rows
  satisfy `last < in_load_release`**, so a `last < in_load_release` guard on the
  UPDATE steps cleanly excludes rows this same release just inserted — with no
  false exclusions.

### Step-by-step mapping (incremental ⇐ PEL)

Run order (helper snapshot first, then inserts, then in-place updates):

1. `load_upi_max_versions_table(dbid)` + `load_max_versions_table()` — snapshot of
   the *original* xref (must precede all mutations).
2. **A1** `incremental_new_versions` ⇐ `populate_pel_tables2`: INSERT a new active
   generation for accessions whose sequence version changed (same upi).
3. **A2/A3** `incremental_new_accessions` ⇐ `populate_pel_tables3` (main + Gap A):
   INSERT active rows for brand-new / previously-fully-deleted accessions and for a
   new sequence variant of an already-active accession.
4. **B1** `incremental_refresh` ⇐ `populate_pel_tables1` (deleted='N' branch):
   in-place UPDATE of matched-unchanged rows (`last`, `taxid`).
5. **B2** `incremental_retire_changed` ⇐ `populate_pel_tables1` (deleted='Y'
   branch): in-place UPDATE retiring matched-but-changed rows.
6. **B3** `incremental_retire_dropped` ⇐ `populate_pel_tables4` (second `UNION ALL`
   arm): in-place UPDATE retiring active rows whose `(ac, upi)` vanished from the
   load. **This is the step with no counterpart in the old Oracle incremental path**
   — the "become inactive" behaviour.

`populate_pel_tables4`'s *first* arm (carry forward all historical `deleted='Y'`
rows) has **deliberately no incremental analogue** — leaving those rows untouched
is the entire performance win.

## Proposed sequence of work (once approved)

1. ✅ **Parity harness** — `tests/xref-incremental-parity/` (`schema.sql`, `seed.sql`,
   `run.sh`). Builds two DBs from one seed, runs FULL in one and INCREMENTAL in the
   other, diffs the resulting xref state on the business key.
2. ✅ **Port + deactivation** — the Oracle `incremental1/2/3` are replaced by
   `incremental_new_versions`, `incremental_new_accessions`, `incremental_refresh`,
   `incremental_retire_changed`, `incremental_retire_dropped` (each the in-place
   analogue of a `populate_pel_tables*` step), orchestrated by the rewritten
   `load_xref_incremental`. The `version_i` divergence is resolved (we follow the
   PEL rule), and case 5 deactivation is covered by `incremental_retire_dropped`.
3. ✅ **Parity validated** on the fixture: full and incremental produce byte-identical
   xref state across cases 1/2/5/D/Gap A + reactivation + deleted-history carry-forward.
   Confirmed incremental leaves unchanged rows physically in place (stable `id`) while
   full reassigns every `id`.
4. ✅ **Selection logic** — new `release.get_load_release_type(dbid)` returns `'F'`
   for a database with no prior release (bootstrap) and `'I'` otherwise.
   `rnc_update.prepare_releases` gains an `'A'` (auto) mode that applies that choice
   per database; an explicit `'F'`/`'I'` still forces a type. `run.run()` calls
   `prepare_releases('A')` by default, and `rnac release run --force-full` forces
   `'F'` for every database (the escape hatch, replacing the reliance on a per-row
   `force_load` flag that nothing set). Tested by:
   - `tests/rnacentral/release/run_test.py` — `run()` issues `prepare_releases('A')`
     by default and `prepare_releases('F')` under `--force-full`.
   - `tests/xref-incremental-parity/selection.sh` — DB-level assertions that
     `get_load_release_type` and `prepare_releases('A'|'F')` assign the right types.
5. ⬜ **Staging validation + benchmark** — run the harness pattern against a real
   staging copy per database; measure wall-clock vs full; watch bloat/autovacuum on
   the mutated live partitions. **Do this before relying on the auto default in
   production**; until then, `--force-full` keeps the old behaviour.
6. ⬜ **Roll out** — drop `--force-full` once staging confirms parity + speedup.

### Validating the harness locally

    PGHOST=127.0.0.1 PGPORT=5433 PGUSER=postgres ./tests/xref-incremental-parity/run.sh

(any Postgres you can create databases in; it creates/drops `parity_full` and
`parity_incr`). Extend `seed.sql` with more edge cases and re-run to widen coverage
before enabling incremental in step 4.

## Out of scope (noted)

The staging pgloader step still truncate-reloads the full CSV because the parsers
emit a full dump. Making *that* incremental would require parsers to diff against
DB state — larger change, separate effort.
