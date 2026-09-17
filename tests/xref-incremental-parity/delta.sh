#!/usr/bin/env bash
#
# Tests DELTA-mode xref loading (docs/incremental-parsing.md): the load contains
# only the new/changed records, deletions come from an explicit load_deletions list,
# and rows that are merely ABSENT from the load must stay active.
#
# Requires a running PostgreSQL you can create databases in (libpq env vars), e.g.
#
#   PGHOST=127.0.0.1 PGPORT=5433 PGUSER=postgres ./delta.sh
#
# It creates/drops a database named delta_test.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
XDIR="$ROOT/database_functions/rnc_load_xref"
IDIR="$ROOT/database_functions/rnc_load_xref_incremental"

P () { psql -v ON_ERROR_STOP=1 -q "$@"; }

DB=delta_test
P -d postgres -c "drop database if exists $DB;"
P -d postgres -c "create database $DB;"
P -d "$DB" -f "$HERE/schema.sql"

# load_xref_delta needs the two version-snapshot helpers plus the incremental steps.
for f in load_upi_max_versions_table load_max_versions_table; do
  P -d "$DB" -c "set search_path=rnacen,public;" -f "$XDIR/$f.sql"
done
for f in "$IDIR"/*.sql; do P -d "$DB" -c "set search_path=rnacen,public;" -f "$f"; done

P -d "$DB" <<'SQL'
set search_path = rnacen, public;
set client_min_messages = warning;

insert into rnacen.rna (urs, len) values
  ('UPI_K',100),('UPI_C',100),('UPI_N',100),('UPI_D',100),('UPI_O',100);

-- Initial state from release 1 (created=last=1).
insert into rnacen.xref (ac, dbid, version, version_i, urs, created, last, deleted, taxid) values
  ('ACC_KEEP', 1, 1, 1, 'UPI_K', 1, 1, 'N', 9606),  -- unchanged AND absent from delta load
  ('ACC_CHG',  1, 1, 1, 'UPI_C', 1, 1, 'N', 9606),  -- version changes
  ('ACC_DEL',  1, 1, 1, 'UPI_D', 1, 1, 'N', 9606),  -- explicitly deleted
  ('ACC_OLD',  1, 1, 1, 'UPI_O', 1, 1, 'Y', 9606);  -- already deleted

-- Delta load: ONLY the changed + new records.
insert into rnacen.load_retro_tmp (in_dbid, in_load_release, in_ac, in_version, in_taxid, comparable_prot_upi) values
  (1, 2, 'ACC_CHG', 2, 9606, 'UPI_C'),
  (1, 2, 'ACC_NEW', 1, 9606, 'UPI_N');

-- Explicit deletion list (from the parser's manifest diff).
insert into rnacen.load_deletions (database, accession) values ('TESTDB', 'ACC_DEL');

select rnc_load_xref_incremental.load_xref_delta(1::bigint, 1::bigint);

DO $$
BEGIN
  -- The whole point: an unchanged row that is absent from the delta load stays active.
  ASSERT (select deleted from xref where ac='ACC_KEEP') = 'N', 'ACC_KEEP must stay active';
  ASSERT (select last    from xref where ac='ACC_KEEP') = 1,   'ACC_KEEP must be untouched (last unchanged)';

  -- Explicit deletion retired ACC_DEL.
  ASSERT (select deleted from xref where ac='ACC_DEL') = 'Y', 'ACC_DEL must be retired';

  -- Changed accession: old version retired, new active generation inserted.
  ASSERT (select deleted from xref where ac='ACC_CHG' and version=1) = 'Y', 'old ACC_CHG must retire';
  ASSERT (select deleted from xref where ac='ACC_CHG' and version=2) = 'N', 'new ACC_CHG must be active';

  -- New accession inserted and active.
  ASSERT (select deleted from xref where ac='ACC_NEW') = 'N', 'ACC_NEW must be active';

  -- Already-deleted history untouched.
  ASSERT (select count(*) from xref where ac='ACC_OLD' and deleted='Y') = 1, 'ACC_OLD unchanged';

  -- Exactly one active row per still-present accession; ACC_KEEP never duplicated.
  ASSERT (select count(*) from xref where deleted='N') = 3, 'expected 3 active rows';
END $$;
SQL

P -d postgres -c "drop database if exists $DB;" >/dev/null
echo ">>> DELTA OK: absent rows stay active; only explicit deletions retire <<<"
