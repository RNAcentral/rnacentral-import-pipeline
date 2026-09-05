#!/usr/bin/env bash
#
# Tests the FULL-vs-INCREMENTAL release-type selection (step 4 of
# improvement/load-only-new-data):
#
#   * release.get_load_release_type(dbid) -> 'F' for a database with no prior
#     release, 'I' for one that has been loaded before.
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
insert into rnacen.rnc_database values (1, 'DBONE'), (2, 'DBTWO');

create table rnacen.rnc_release (
  id bigint primary key, dbid smallint, release_date date, release_type char(1),
  status char(1), "timestamp" timestamp default now(), userstamp text, descr text, force_load char(1)
);
-- DBONE has been loaded before (a completed release); DBTWO never has.
insert into rnacen.rnc_release (id, dbid, release_type, status) values (1, 1, 'F', 'D');

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
  ASSERT release.get_load_release_type(1) = 'I', 'db with prior release should be INCREMENTAL';
  ASSERT release.get_load_release_type(2) = 'F', 'db with no prior release should be FULL';
  ASSERT release.get_load_release_type(3) = 'F', 'unknown db should default to FULL';
END $$;

-- 2. prepare_releases('A') picks per database.
select rnc_update.prepare_releases('A');
DO $$
BEGIN
  ASSERT (select release_type from rnc_release where dbid = 1 and status = 'L') = 'I',
         'auto mode should give DBONE an INCREMENTAL release';
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
--    instead of INCREMENTAL; a database with no manifest still selects INCREMENTAL.
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
SQL

P -d postgres -c "drop database if exists $DB;" >/dev/null
echo ">>> SELECTION OK: get_load_release_type / prepare_releases behave as specified <<<"
