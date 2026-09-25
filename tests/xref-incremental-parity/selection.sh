#!/usr/bin/env bash
#
# Tests the FULL-vs-DELTA release-type selection (step 4 of
# improvement/load-only-new-data):
#
#   * release.get_load_release_type(dbid) -> 'F' with no prior release or no
#     import manifest, 'D' once a manifest exists for the database -- except
#     HGNC, paused on 'F' regardless of its manifest (see the function's comment),
#     and WormBase, which follows ENA's manifest.
#   * rnc_update.prepare_releases('A') applies that choice per database.
#   * rnc_update.prepare_releases('F') forces FULL for every database.
#
# Requires a running PostgreSQL you can create databases in (libpq env vars), e.g.
#
#   PGHOST=127.0.0.1 PGPORT=5433 PGUSER=postgres ./selection.sh
#
# It creates/drops a database named selection_test.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"

P () { psql -v ON_ERROR_STOP=1 -q "$@"; }

DB=selection_test
P -d postgres -c "drop database if exists $DB;"
P -d postgres -c "create database $DB;"

P -d "$DB" <<'SQL'
set client_min_messages = warning;
create schema rnacen;
create schema rnc_update;
create schema release;
set search_path = rnacen, public;

create table rnacen.rnc_database (id smallint primary key, descr text);
insert into rnacen.rnc_database values (1, 'DBONE'), (2, 'DBTWO'), (4, 'HGNC'), (5, 'WORMBASE'), (6, 'ENA');

create table rnacen.rnc_release (
  id bigint primary key, dbid smallint, release_date date, release_type char(1),
  status char(1), "timestamp" timestamp default now(), userstamp text, descr text, force_load char(1)
);
-- DBONE and HGNC have both been loaded before (a completed release); DBTWO never has.
insert into rnacen.rnc_release (id, dbid, release_type, status) values (1, 1, 'F', 'D'), (2, 4, 'F', 'D'), (3, 5, 'F', 'D');

create table rnacen.load_rnacentral_all (database varchar(40));
insert into rnacen.load_rnacentral_all values ('DBONE'), ('DBTWO');
SQL

# Load the functions under test straight from the source tree.
P -d "$DB" -c "set search_path=rnacen,public;" -f "$ROOT/database_functions/release/get_load_release_type.sql"
P -d "$DB" -c "set search_path=rnacen,public;" -f "$ROOT/database_functions/rnc_update/create_release.sql"
P -d "$DB" -c "set search_path=rnacen,public;" -f "$ROOT/database_functions/rnc_update/prepare_releases.sql"

P -d "$DB" <<'SQL'
set search_path = rnacen, public;
set client_min_messages = warning;

-- 1. get_load_release_type reflects history.
DO $$
BEGIN
  ASSERT release.get_load_release_type(1) = 'F', 'db with prior release but no manifest should be FULL';
  ASSERT release.get_load_release_type(2) = 'F', 'db with no prior release should be FULL';
  ASSERT release.get_load_release_type(3) = 'F', 'unknown db should default to FULL';
  ASSERT release.get_load_release_type(4) = 'F', 'HGNC should be FULL';
END $$;

-- 2. prepare_releases('A') picks per database.
select rnc_update.prepare_releases('A');
DO $$
BEGIN
  ASSERT (select release_type from rnc_release where dbid = 1 and status = 'L') = 'F',
         'auto mode should give DBONE a FULL release (no manifest)';
  ASSERT (select release_type from rnc_release where dbid = 2 and status = 'L') = 'F',
         'auto mode should give DBTWO a FULL release';
END $$;

-- 3. prepare_releases('F') forces FULL for all (escape hatch). Clear the pending
--    releases first, otherwise prepare_releases short-circuits on status='L'.
delete from rnc_release where status = 'L';
select rnc_update.prepare_releases('F');
DO $$
BEGIN
  ASSERT (select count(*) from rnc_release where status = 'L' and release_type <> 'F') = 0,
         'forced mode should give every database a FULL release';
  ASSERT (select count(*) from rnc_release where status = 'L') = 2,
         'forced mode should create one release per loaded database';
END $$;

-- 4. A delta-parsed database (one that has an import manifest) selects DELTA
--    instead of FULL; a database with no manifest still selects FULL.
create table rnacen.pipeline_tracking_import (
  database text, accession text, signature text,
  updated_at timestamptz default now(), primary key (database, accession)
) partition by list (database);
create table rnacen.pipeline_tracking_import_dbone
  partition of rnacen.pipeline_tracking_import for values in ('DBONE');
insert into rnacen.pipeline_tracking_import (database, accession, signature)
  values ('DBONE', 'x', 'sig');
DO $$
BEGIN
  ASSERT release.get_load_release_type(1) = 'D', 'db with a manifest selects DELTA';
END $$;

-- 5. HGNC is the exception: its manifest is real (its parser was the first one
--    built this way), but auto-DELTA for it is paused, so it stays FULL even
--    with a manifest present.
create table rnacen.pipeline_tracking_import_hgnc
  partition of rnacen.pipeline_tracking_import for values in ('HGNC');
insert into rnacen.pipeline_tracking_import (database, accession, signature)
  values ('HGNC', 'x', 'sig');
DO $$
BEGIN
  ASSERT release.get_load_release_type(4) = 'F', 'HGNC stays FULL even with a manifest';
END $$;

-- 6. WormBase has no manifest of its own: its xrefs come out of the ENA parse, so
--    it follows ENA's. FULL would retire every WormBase record in the ENA sources
--    the delta skipped.
DO $$
BEGIN
  ASSERT release.get_load_release_type(5) = 'F', 'WormBase is FULL while ENA has no manifest';
END $$;
create table rnacen.pipeline_tracking_import_ena
  partition of rnacen.pipeline_tracking_import for values in ('ENA');
insert into rnacen.pipeline_tracking_import (database, accession, signature)
  values ('ENA', 'x', 'sig');
DO $$
BEGIN
  ASSERT release.get_load_release_type(5) = 'D', 'WormBase follows ENA into DELTA';
END $$;

-- 7. Every ENA run releases WormBase, even with no WormBase rows staged: a run can
--    delete WormBase records without loading any.
delete from rnc_release where status = 'L';
delete from load_rnacentral_all;
insert into load_rnacentral_all values ('ENA');
select rnc_update.prepare_releases('A');
DO $$
BEGIN
  ASSERT (select release_type from rnc_release where dbid = 5 and status = 'L') = 'D',
         'an ENA run must create a DELTA WormBase release';
END $$;
SQL

# ...and release run must execute it: TO_RELEASE is read straight from run.py.
TO_RELEASE="$(cd "$ROOT" && python3 -c 'from rnacentral_pipeline.rnacentral.release import run; print(run.TO_RELEASE)')"
released="$(P -d "$DB" -At -c "set search_path=rnacen,public;" -c "$TO_RELEASE" | cut -d'|' -f1 | tr '\n' ' ')"
if [[ " $released " != *" 5 "* ]]; then
  echo "TO_RELEASE must pick up WormBase's release on an ENA run (got: $released)" >&2
  exit 1
fi

P -d postgres -c "drop database if exists $DB;" >/dev/null
echo ">>> SELECTION OK: get_load_release_type / prepare_releases behave as specified <<<"
