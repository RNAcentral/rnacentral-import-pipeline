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

## Workflow

```
fetch_directory ─▶ .ncr chunks
     │
     ├─▶ ena_signatures(chunk)         (cheap, parallel: record.id,signature)
     │        └─ collect ─▶ ena_delta_diff(db)   (single reduce, DB-side)
     │                          ├─▶ to_parse.txt   (accessions: new + changed)
     │                          ├─▶ deletions.csv  (ENA,accession for dropped)
     │                          └─▶ manifest.csv   (ENA,accession,signature: all current)
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
`rnac manifest apply`. No import-data.nf change is needed — ENA just emits the same
two side-channel files HGNC does, keyed by the `database` column.

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

## DB-side diff

`ena_delta_diff` `COPY`s the collected `(accession, signature)` into a temp table
and computes, in SQL against `pipeline_tracking_import` scoped to ENA:

- **to_parse** — temp rows with no matching stored row, or a differing signature;
- **deletions** — stored ENA rows absent from temp;
- **manifest.csv** — every temp row (the full new manifest).

This keeps the millions-row set difference in Postgres instead of Python.

## Load side — unchanged, already generic

Nothing new is needed on the load side. `get_load_release_type(dbid)` already
returns `D` for *any* database that has manifest rows, and `load_xref_delta` +
`incremental_retire_explicit` + `load_deletions` already retire only explicitly
listed accessions. ENA becomes a delta database purely by producing a manifest.
`rnac release run --force-full` remains the escape hatch (forces `F`, ignores the
manifest) and reconciles any signature/mapping drift.

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
  without its signature changing. `--force-full` periodically reconciles.
- **Not yet benchmarked at production volume.** The real question — what fraction of
  ENA changes between snapshots — is unmeasured; that fraction is the whole payoff.
