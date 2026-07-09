# Database functions

This directory is the **source of truth** for the in-database PL/pgSQL functions and
procedures the release pipeline relies on. One file per function:

```
database_functions/<schema>/<function_name>.sql
```

Each file is a complete `CREATE OR REPLACE FUNCTION ...` (as produced by
`pg_get_functiondef`). Filenames are the bare function name — the directory supplies the
schema. Overloaded functions (same name, different arg types) are **not** supported by this
layout; `fetch_functions.sh` aborts if it finds any.

## How to change a function

1. Edit the relevant `database_functions/<schema>/<func>.sql` file.
2. Commit it.
3. It deploys on the next release run: `run()` calls `functions.apply()`, which applies
   every file whose contents changed since it was last applied. To deploy ad hoc:

   ```
   rnac release functions            # apply changed functions
   rnac release functions --dry-run  # show what would change, apply nothing
   ```

Do **not** hand-edit functions directly in the database — that creates drift (see below).
There are no more inline `CREATE OR REPLACE` patches in `run.py`; that mechanism was
retired in favour of this directory + the applier.

## How deployment works

- `rnacentral_pipeline/rnacentral/release/functions.py` discovers the files, sha256s each,
  and in a **single transaction** applies only those whose hash differs from what was last
  applied. It records `(schema_name, function_name, sha256, applied_at)` in
  `rnacen.applied_functions`.
- `CREATE OR REPLACE` for PL/pgSQL is order-independent (bodies aren't resolved against
  other objects until run time), so files apply in plain `(schema, name)` order — no
  dependency graph.
- **`rnc_test` is tracked here but never deployed** (test fixtures), via
  `DEPLOY_EXCLUDE_SCHEMAS` in `functions.py`. Add a schema to that set to track-but-not-deploy.

## Refreshing from the live database / drift check

`fetch_functions.sh` re-dumps functions from a database into this exact layout:

```
PGDATABASE=<conninfo> ./fetch_functions.sh database_functions \
    database release rnacen rnc_load_rna rnc_load_xref rnc_load_xref_incremental \
    rnc_logging rnc_test rnc_update stats expunge trembl_utils upi
```

Then `git diff` shows any drift between the repo and the live database. (`public` is
intentionally omitted — its functions are extension-owned.)

## Index DDL is *not* here

Index creation used in the release (e.g. `load_md5_new_sequences$in_md5`) lives in
`run.py`, not this directory — the applier manages functions only.
