# Incremental (delta) parsing — ENA

**Branch:** `improvement/load-only-new-data`
**Status:** parser-side implementation. Builds on the generic manifest machinery
([manifest.py](../rnacentral_pipeline/databases/manifest.py),
[incremental-parsing.md](incremental-parsing.md)) first used for HGNC.

## Why ENA is not a copy of HGNC

HGNC is a single ~8,400-record JSON file parsed in one process, so the manifest
approach there is: load the whole file into a dict, hash each record, diff against
the stored manifest, map only the changed handful. None of that shape survives for
ENA:

| | HGNC | ENA |
|---|---|---|
| Input | one small JSON | millions of EMBL records |
| Fetch | one `wget` | `rsync` of an NFS snapshot across `con/std/wgs/tls/tsa` |
| Processing | one process | split into 50k-record `.ncr` chunks, **many parallel** `process_file` tasks |
| Expensive step | per-record HTTP/DB in mapping | `ribotyper` (rRNA typing) + parse, per sequence |
| Manifest size | ~8,400 rows | millions of rows |

Three consequences drive the design below:

1. **No single process sees all records.** `dropped = old − new` is a global set
   difference, but ENA records are spread across parallel chunks. The diff must be
   a single reduce over every chunk's signatures, not a per-process step.
2. **The manifest is huge.** A per-row `INSERT` loop (fine for HGNC) is hopeless
   for millions of rows, and holding both the old and new signature maps in Python
   to diff them is memory-hungry. Both the diff and the manifest upsert are pushed
   into Postgres.
3. **The win is `ribotyper` + parse, not the fetch.** `rsync` and the signature
   pass still read every record. The saving is skipping `ribotyper` and the mapping
   for records that did not change — worthwhile only because a chunk with no changed
   records is dropped entirely (no `ribotyper` run at all).

## The record key

`rnac ena parse` iterates with `SeqIO.parse(path, "embl")` and uses `record.id`
as the accession ([ena/helpers.py](../rnacentral_pipeline/databases/ena/helpers.py)
`accession()` → `record.id`). BioPython sets `record.id` to
`primary_accession.sequence_version` (e.g. `AB123456.2`). That is exactly the
`Entry.accession` the loader stores as `xref.ac`, and exactly what the explicit
deletion path (`incremental_retire_explicit`, `d.accession = u.ac`) matches on. So
the manifest, the to-parse filter, and the deletion list all key on `record.id`.

Because the version is part of the key, an ENA sequence-version bump appears as a
new accession (`.N`) plus a dropped one (`.N-1`) — the correct behaviour when xref
identity includes the version.

To guarantee the signature/filter keys are byte-identical to the parser's, the
signature and filter passes derive the key with the *same* `SeqIO` record.id
logic, per raw record block, rather than re-implementing accession parsing.

## Signature

`signature = sha256(raw EMBL record text)` — the whole `ID … //` block. Hashing the
raw text can never produce a false "unchanged" (the dangerous direction); at worst a
cosmetic reformat causes a needless re-parse (safe, just slower). This matches the
HGNC rule "whole record, so any change is caught".

## Source signature — the prefilter

Records are the finest level, but most of the work happens before we ever see one:
rsyncing the archives, decompressing them, concatenating and splitting them, and the
signature pass itself. So there is a coarser level above the records:

`source signature = sha256(sorted "relative path, size, mtime" of every archive under
the source)`, where a **source** is one project subdirectory of `wgs/`, `tls/` or
`tsa/`, or the whole of `con/` or `std/`. It is computed from the directory listing
alone, so it costs a `find` and no data is read; a source whose signature has not
moved is not copied, not decompressed, not split and not signatured.

The two parts of the tree behave completely differently, which is why both levels
earn their place:

- **`wgs/`, `tls/`, `tsa/`** are per-project archives that are written once and then
  left alone — `wgs/public/aaa/AAAABA02.ncr.gz` was last modified in October 2024.
  Whole subdirectories are identical between snapshots, and the source signature
  skips them outright.
- **`con/` and `std/`** are numbered slices of a regenerated dump (`STD_XXX_42`), so
  the slicing shifts whenever anything upstream changes. Of the 120 `.ncr.gz`
  archives present under both published monthly snapshots (`snapshot_20260726` and
  `snapshot_20260826`), **none** kept the same size, let alone the same bytes; July
  had 333 of them, August 120. The source signature never skips these, and only the
  per-record key finds anything: hashing records over that same snapshot pair finds
  96-97% of them byte-identical (`STD_MAM_1`: 18832 of 19488; `STD_ENV_1`: 15908 of
  16337).

Hashing our own chunks would work at neither level — `fetch_directory` splits each
source with `split-ena`, so a chunk boundary moves whenever a record ahead of it is
added or removed.

### Skipping a source without losing its records

The record diff retires accessions by absence from the new signature set, so a
skipped source — which contributes no signatures — would retire everything it holds.
Two things prevent that:

* `pipeline_tracking_import.source_id` records which source each accession came from,
  pointing at `pipeline_tracking_import_files(id, database, path, signature)`. It is a
  surrogate id rather than the path itself because at ENA row counts 8 bytes per row
  beats repeating a ~60 byte path in every one of them.
* `ena source-diff` writes `scanned.txt` — the sources fetched this run plus the ones
  that have gone from the snapshot — and `diff_via_polars` will only delete records
  belonging to those. A skipped source is in neither set, so its records are left
  exactly as they are.

This is also why `--force_full_import` forces every source to be fetched: that run
releases with `F`, which retires by absence from the load.

Chunks carry their origin in their own name: `fetch_directory` handles one source at a
time and names its output `<label>-chunk<N>.ncr`, where the label is the first 16 hex
digits of `sha1(source path)`. That is the only thing tying a chunk back to the
directory it came from once the fetch is over.

## Workflow

```
list_subdirs ─▶ stat_sources (parallel: path,signature from the listing only)
     └─ collect ─▶ ena_source_diff(db)
                        ├─▶ to_fetch.txt  (sources new or changed: the only fetches)
                        ├─▶ scanned.txt   (those plus vanished: retirable sources)
                        └─▶ sources.csv   (ENA,path,signature: all current)

fetch_directory(to_fetch) ─▶ .ncr chunks, one source at a time
     │
     ├─▶ ena_signatures(chunk)         (cheap, parallel: record.id,signature)
     │        └─ collect ─▶ ena_delta_diff(db, scanned.txt)   (single reduce)
     │                          ├─▶ to_parse.txt   (accessions: new + changed)
     │                          ├─▶ deletions.csv  (ENA,accession for dropped)
     │                          └─▶ manifest.csv   (ENA,accession,signature,source)
     │
     └─▶ process_file(chunk, to_parse, metadata)   (parallel)
              rnac ena filter --only to_parse.txt chunk → filtered.ncr
              (filtered.ncr empty ⇒ emit nothing, no ribotyper)
              ena2fasta ▸ ribotyper ▸ rnac ena parse filtered.ncr
     emit: data csvs + manifest.csv + deletions.csv
```

`manifest.csv` and `deletions.csv` ride the existing generic wiring
([import-data.nf](../import-data.nf)): `deletions.csv` is loaded into
`load_deletions` via `deletions.ctl`; `manifest.csv` is promoted post-release by
`rnac manifest apply`. `sources.csv` is promoted the same way and for the same
reason, by `rnac manifest apply-sources` — a run that never released must not leave
the files table claiming those sources are up to date.

### Bootstrap (first delta run)

ENA already has release history, so before its first manifest exists
`get_load_release_type` returns `I` (incremental, absence = deleted) and the diff
sees an empty prior manifest. Rather than materialise a to-parse set containing
*every* accession (millions, loaded into every chunk's filter), `ena_delta_diff`
writes the sentinel `__ALL__` as the sole line of `to_parse.txt` when there is no
prior ENA manifest. `rnac ena filter` treats `__ALL__` as "keep everything" and
copies the chunk through unchanged. Deletions are empty on bootstrap. The run parses
everything (as today) and seeds the manifest; the next run sees the manifest, loads
in `D` mode, and only changed records are parsed.

## The diff — on the cluster, in polars

The database is only ever read from. `ena_delta_diff` `COPY`s ENA's stored manifest
out to `stored-manifest.csv` (`manifest.dump_signatures`, one statement against one
partition) and does every join here in polars (`manifest.diff_via_polars`):

- **to_parse** — new rows with no stored match, or a differing signature;
- **deletions** — stored ENA rows absent from the new set, restricted to rows whose
  source was scanned this run;
- **manifest.csv** — every new row, with the source it was read from.

Both joins run under polars' streaming engine, so the millions-row set difference
never has to fit in memory, and a diff that lands while the database is busy cannot
contend with anything: no temp tables, no server-side sorting, no writes.

## Load side — unchanged, already generic

Nothing new is needed on the load side. `get_load_release_type(dbid)` already
returns `D` for *any* database that has manifest rows, and `load_xref_delta` +
`incremental_retire_explicit` + `load_deletions` already retire only explicitly
listed accessions. ENA becomes a delta database purely by producing a manifest.
`nextflow run main.nf --force_full_import` remains the escape hatch: the diff is
skipped, every chunk is parsed, and the release forces `F`, so a run reconciles any
signature/mapping drift.

## Safety / rollout

Until an ENA manifest exists the system is unchanged: full parse every time + `I`
load (absence = deleted). A partial or failed rollout cannot corrupt data because
`store_signatures`/`apply` run only after the release commits, and the DELTA load
path is selected only when a manifest is actually present.

## Known limitations / follow-ups

- **Signature pass cost.** v1 derives `record.id` via `SeqIO` per raw block for
  exactness; it is CPU-only and parallel per chunk but re-reads every record. A Rust
  pass (mirroring `utils/split-ena`) is the obvious optimisation once correctness is
  confirmed on real data.
- **Mapping drift.** Like HGNC, the signature hashes only the raw record, so a record
  whose *mapping* depends on RNAcentral's own data could in principle map differently
  without its signature changing. `--force_full_import` periodically reconciles.
- **Not yet benchmarked at production volume.** The real question — what fraction of
  ENA changes between snapshots — is unmeasured; that fraction is the whole payoff.
