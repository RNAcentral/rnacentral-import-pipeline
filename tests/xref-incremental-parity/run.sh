#!/usr/bin/env bash
#
# Parity harness for the incremental xref load (improvement/load-only-new-data).
#
# Builds two throwaway databases from the same seed, runs the FULL Partition
# Exchange Load in one and the INCREMENTAL load in the other, then diffs the
# resulting xref state. They must be byte-identical on the business key
# (id is intentionally excluded: a full release reassigns every id, an
# incremental one preserves them).
#
# Requires a running PostgreSQL you can create databases in. Point it at one via
# the standard libpq env vars (PGHOST/PGPORT/PGUSER/PGPASSWORD), e.g.
#
#   PGHOST=127.0.0.1 PGPORT=5433 PGUSER=postgres ./run.sh
#
# It creates/drops databases named parity_full and parity_incr.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
XDIR="$ROOT/database_functions/rnc_load_xref"
IDIR="$ROOT/database_functions/rnc_load_xref_incremental"

P () { psql -v ON_ERROR_STOP=1 -q "$@"; }

# Full-path functions the PEL load depends on (load order is irrelevant:
# CREATE OR REPLACE for plpgsql resolves bodies at run time).
FULL_FNS=(
  load_upi_max_versions_table load_max_versions_table prepare_pel_tables
  populate_pel_tables1 populate_pel_tables2 populate_pel_tables3 populate_pel_tables4
  do_pel_exchange load_xref do_checks
)

setup_db () {
  local db="$1"
  P -d postgres -c "drop database if exists $db;"
  P -d postgres -c "create database $db;"
  P -d "$db" -f "$HERE/schema.sql"
  for f in "${FULL_FNS[@]}"; do P -d "$db" -c "set search_path=rnacen,public;" -f "$XDIR/$f.sql"; done
  for f in "$IDIR"/*.sql;   do P -d "$db" -c "set search_path=rnacen,public;" -f "$f"; done
  P -d "$db" -f "$HERE/seed.sql"
}

DUMP="set search_path=rnacen,public;
      select ac,dbid,version,version_i,urs,created,last,deleted,coalesce(taxid,-1) as taxid
      from xref order by ac,urs,version,deleted,created;"

setup_db parity_full
setup_db parity_incr

P -d parity_full -c "set search_path=rnacen,public; select rnc_load_xref.load_xref(1::bigint, 1::bigint);" >/dev/null
P -d parity_incr -c "set search_path=rnacen,public; select rnc_load_xref_incremental.load_xref_incremental(1::bigint, 1::bigint);" >/dev/null

full="$(mktemp)"; incr="$(mktemp)"
P -d parity_full -At -F$'\t' -c "$DUMP" > "$full"
P -d parity_incr -At -F$'\t' -c "$DUMP" > "$incr"

if diff -u "$full" "$incr"; then
  echo ">>> PARITY OK: full and incremental produced identical xref state <<<"
else
  echo ">>> PARITY FAILED: see diff above <<<"; exit 1
fi
