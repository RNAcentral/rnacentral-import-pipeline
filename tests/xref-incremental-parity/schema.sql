-- Minimal fixture reproducing the RNAcentral xref loading structures.
set client_min_messages = warning;

create schema if not exists rnacen;
create schema if not exists rnc_load_xref;
create schema if not exists rnc_load_xref_incremental;
set search_path = rnacen, public;

-- reference tables (FK targets used by do_pel_exchange)
create table rnacen.rnc_database (id smallint primary key, descr text);
insert into rnacen.rnc_database values (1, 'TESTDB');

create table rnacen.rnc_release (
  id bigint primary key, dbid smallint, release_date date, release_type char(1),
  status char(1), "timestamp" timestamp default now(), userstamp text, descr text, force_load char(1)
);
insert into rnacen.rnc_release (id, dbid, release_type, status) values
  (1, 1, 'F', 'D'),   -- previous release (done)
  (2, 1, 'F', 'L');   -- current release being loaded

create table rnacen.rna (upi varchar(30) primary key, len int);

-- xref: declarative partitioning dbid LIST -> deleted LIST
create table rnacen.xref (
  id bigserial,
  ac varchar(200), dbid smallint, version bigint, version_i bigint, upi varchar(30),
  created bigint, last bigint, deleted varchar(1), taxid bigint,
  "timestamp" timestamp default clock_timestamp(), userstamp text default current_user
) partition by list (dbid);

create table rnacen.xref_p1 partition of rnacen.xref for values in (1) partition by list (deleted);
create table rnacen.xref_p1_not_deleted partition of rnacen.xref_p1 for values in ('N');
create table rnacen.xref_p1_deleted     partition of rnacen.xref_p1 for values in ('Y');
create unique index "xref_p1_not_deleted$id" on rnacen.xref_p1_not_deleted (id);
create unique index "xref_p1_deleted$id"     on rnacen.xref_p1_deleted (id);

-- load tables (subset used by the xref path)
create unlogged table rnacen.load_retro_tmp (
  in_dbid smallint, in_load_release integer, in_crc64 varchar(16), in_len integer,
  in_seq_short varchar(4000), in_seq_long text, in_ac varchar(200), in_version bigint,
  in_md5 varchar(32), in_taxid integer, comparable_prot_upi varchar(30)
);
create unlogged table rnacen.load_upi_max_versions (ac varchar(200), dbid smallint, max_version_i bigint, upi varchar(30));
create unlogged table rnacen.load_max_versions (ac varchar(200), dbid smallint, max_version_i bigint);
