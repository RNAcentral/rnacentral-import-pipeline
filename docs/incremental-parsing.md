# Incremental (delta) parsing — design

**Branch:** `improvement/load-only-new-data`
**Status:** design + parser-side foundation. First target: **HGNC**. Load-side
delta mode still to build.

## Problem

Parsing re-processes a database's entire input every import, even when almost
nothing changed. For **HGNC** an import took ~45 min to parse while only 8 records
had actually changed. The cost is per-record work in the mapping step
([hgnc/parser.py](../rnacentral_pipeline/databases/hgnc/parser.py) →
[helpers.py](../rnacentral_pipeline/databases/hgnc/helpers.py)):

- `ensembl_sequence` — a live HTTP GET to `rest.ensembl.org` per Ensembl-mapped record;
- `refseq_id_to_urs`, `gtrnadb_to_urs`, `md5_to_urs`, `urs_to_sequence` — per-record DB queries.

All redundant when the record is unchanged since the last import.

## Idea

Identify what changed *before* the expensive parse, and only parse (and load) the
changed records; mark vanished records inactive explicitly.

1. **Cheap manifest pass** — read the raw input, emit `accession -> signature`
   where signature is a hash of the raw record. No DB, no HTTP.
2. **Diff** against the manifest stored from the last successful import →
   `NEW`, `CHANGED`, `DROPPED`, `UNCHANGED`.
3. **Full parse only `NEW`+`CHANGED`** (the 8 records) — the only ones that hit
   Ensembl / the DB.
4. **`DROPPED`** → an explicit deletion list.
5. **Persist the new manifest** — but only once the load for the database has
   succeeded (see transactionality).

First run (no manifest) → everything is `NEW` → a full parse that also seeds the
manifest. Same bootstrap-then-incremental shape as the loading side.

## The load-side consequence (important)

The incremental loader committed on this branch marks a sequence inactive **by its
absence** from the load (`incremental_retire_dropped`, mirroring
`populate_pel_tables4`). A delta parse emits only the *changed* subset, so that
absence rule would wrongly retire every record we didn't re-emit.

Delta loading therefore needs a **third release mode** alongside FULL (`F`) and
INCREMENTAL (`I`):

- **DELTA (`D`)** — upsert the changed set with the existing incremental steps
  (`incremental_new_versions`, `incremental_new_accessions`, `incremental_refresh`,
  `incremental_retire_changed` — all of which only act on rows *present* in the
  load), and replace `incremental_retire_dropped` with a step that retires only the
  accessions in an explicit `load_deletions` staging table.

So `load_xref_delta` reuses 4 of the 5 incremental functions and swaps the one
absence-based step for an explicit-deletion one.

### "Absent rows are never deleted" is a database-wide switch

`incremental_retire_dropped` (and `populate_pel_tables4` in the full path) is the
*only* place a row is retired for being absent. Gating that single step per database
fully guarantees "absent ≠ deleted" for that database — nothing else deletes by
absence. The switch is self-managing: **a database has a manifest iff it is
delta-parsed**, so `release.get_load_release_type(dbid)` becomes

- no prior release → `F` (bootstrap full load),
- a manifest exists for the database → `D` (delta: deletions only from the explicit
  list; absent rows are left active),
- otherwise → `I` (legacy full dump: absence = deleted).

This **cannot** be flipped globally for every database yet: a legacy database still
ships a full dump each import and relies on absence to detect drops, so turning
absence-retire off for it would mean dropped records are never retired. It is a
per-database switch (default: keep absence-retire; delta databases: off). Once every
database is delta-parsed, absence-retire can be removed entirely.

## Signature & the correctness caveat

Signature = sha256 of the canonicalised raw record (whole record, so any change is
caught). The mapping HGNC produces also depends on RNAcentral's *own* data (e.g.
the NONCODE/RefSeq → URS lookups), so an unchanged raw record could in principle
map differently after the DB changes. Signature-by-raw-hash assumes the mapping is
a pure function of the raw record; a periodic full re-import (the same
`--force-full` escape hatch as the loader) reconciles any drift.

## Manifest storage & transactionality

A single generic table serves every database:

```
rnc_import_manifest(database text, accession text, signature text,
                    updated_at timestamptz, primary key (database, accession))
```

The manifest represents "signatures of what is currently loaded". It must be
updated **only after the database's load/release commits**, so a failed load
leaves it untouched and the next run simply re-emits the same changed set
(idempotent). Concretely the parse step emits, per database:

- the normal CSVs, for `NEW`+`CHANGED` only;
- `deletions.csv` — the `DROPPED` accessions;
- `manifest.csv` — the full new `accession -> signature` map;

and the load step, on success, applies `deletions.csv` and replaces the stored
manifest from `manifest.csv`.

## Pieces

- ✅ **`rnacentral_pipeline/databases/manifest.py`** — signature, diff, manifest
  table read/write, and `write_artifacts`/`apply_artifacts`. Reusable across
  databases. Unit-tested (`tests/databases/manifest_test.py`).
- ✅ **HGNC parser** — `parser.parse(path, db_url, previous_signatures)` returns a
  `ParseResult(entries, signatures, deletions)`, mapping only `NEW`+`CHANGED` (never
  connecting when nothing changed); `rnac hgnc map` reads the manifest and emits
  `manifest.csv` + `deletions.csv`. Unit-tested (`tests/databases/hgnc/delta_test.py`).
- ✅ **Load-side DELTA mode** — `load_xref_delta` (reuses 4 of 5 incremental steps,
  swaps in `incremental_retire_explicit`), `load_deletions` staging + `deletions.ctl`,
  and the `D` branch in `load_release` / `get_load_release_type`. Harness-tested
  (`tests/xref-incremental-parity/delta.sh`, `selection.sh`).
- ✅ **Manifest application** — `rnac manifest apply` upserts the manifest and applies
  deletions; wired in `import-data.nf` to run only after the release commits.
- ✅ **Workflow wiring** — `manifest.csv` is branched out and applied post-release;
  `deletions.csv` rides the normal CSV stream into `load_deletions` via its ctl.
  `hgnc.nf` is unchanged (already emits `*.csv` and has DB access).

### How to test end-to-end on staging

1. First run of HGNC on a DB with no manifest → bootstrap: parses everything, loads
   FULL, writes the manifest. (Same as today, plus the manifest.)
2. Second run with only a few records changed → the parse should map only those
   records (watch the `HGNC delta: N new, M changed, ...` log line), the load runs as
   `D`, and untouched HGNC xrefs stay active. Confirm with the row counts / timings.
3. `--force-full` (on `rnac release run`) still forces a full rebuild if needed.

Until the manifest exists the system safely degrades to full-parse-every-time +
incremental load (correct, just not sped up), so partial wiring cannot corrupt data.

## Generalisation

Everything except the "read the raw records and key them by accession" step is
database-agnostic (`manifest.py`, the DELTA load mode, the workflow plumbing). Each
new database only needs a cheap `raw -> {accession: raw_record}` reader, which most
already have. Roll out per database, slowest first, once HGNC proves the pattern.
```
