#!/usr/bin/env bash
set -euo pipefail

# Usage: ./fetch_functions.sh OUTDIR schema1 schema2 ...
# Dumps every function/procedure in the named schemas to OUTDIR/<schema>/<name>.sql
# (one file per function, full CREATE OR REPLACE). This layout is the committed
# source of truth; re-running over the same OUTDIR lets you `git diff` for drift.
#
# Connection via libpq env vars (PGHOST/PGPORT/PGUSER/PGDATABASE/PGPASSWORD),
# or set CONN to psql flags, e.g. CONN="-h host -U andrew -d rnacen"
CONN="${CONN:-}"

OUTDIR="${1:?need an output dir}"; shift
SCHEMAS=("$@")
[ ${#SCHEMAS[@]} -gt 0 ] || { echo "give at least one schema"; exit 1; }

inlist=$(printf "'%s'," "${SCHEMAS[@]}"); inlist="${inlist%,}"
mkdir -p "$OUTDIR"

list_sql="
SELECT n.nspname, p.proname, p.oid,
       count(*) OVER (PARTITION BY n.nspname, p.proname) AS overloads
FROM pg_proc p
JOIN pg_namespace n ON n.oid = p.pronamespace
WHERE n.nspname IN ($inlist)
  AND p.prokind IN ('f','p')        -- functions + procedures; skips aggregates/window fns
  -- AND NOT EXISTS (SELECT 1 FROM pg_depend d WHERE d.objid=p.oid AND d.deptype='e')  -- skip extension fns
ORDER BY 1,2,3;"

collision=0
while IFS=$'\t' read -r schema name oid overloads; do
  [ -z "$schema" ] && continue
  # Filenames are schema-qualified by directory; overloaded functions (same name,
  # different args) would collide once the OID is dropped. Bail loudly if so.
  if [ "$overloads" -gt 1 ]; then
    echo "ERROR: $schema.$name is overloaded ($overloads signatures); filename would collide" >&2
    collision=1
    continue
  fi
  mkdir -p "$OUTDIR/$schema"
  file="$OUTDIR/$schema/${name}.sql"
  psql $PGDATABASE -Aqt -c "SELECT pg_get_functiondef($oid);" > "$file"
  echo "wrote $file"
done < <(psql $PGDATABASE -Aqt -F$'\t' -c "$list_sql")

[ "$collision" -eq 0 ] || { echo "aborting: overloaded functions need a naming scheme" >&2; exit 1; }
