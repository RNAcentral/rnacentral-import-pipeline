CREATE OR REPLACE FUNCTION rnc_load_xref.prepare_pel_tables()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    --------------------------------------------------------
    -- Recreate deleted table to inherit latest xref 
    --------------------------------------------------------
    execute 'drop table if exists xref_pel_deleted';
    -- execute 'create table if not exists xref_pel_deleted (like xref including defaults including constraints including indexes including comments)';
    execute 'create table if not exists xref_pel_deleted (like xref including defaults including constraints including comments including generated)';

    --------------------------------------------------------
    --  Drop indexes on deleted table
    --------------------------------------------------------
    execute 'drop index if exists xref_pel_deleted$ac';
    execute 'drop index if exists xref_pel_deleted$created';
    execute 'drop index if exists xref_pel_deleted$upi';
    execute 'drop index if exists xref_pel_deleted$taxid';

    --------------------------------------------------------
    --  Drop FK constraints on deleted table
    --------------------------------------------------------
    execute 'alter table xref_pel_deleted drop constraint if exists xref_pel_deleted_fk1';
    execute 'alter table xref_pel_deleted drop constraint if exists xref_pel_deleted_fk2';
    execute 'alter table xref_pel_deleted drop constraint if exists xref_pel_deleted_fk3';
    execute 'alter table xref_pel_deleted drop constraint if exists xref_pel_deleted_fk4';

    --------------------------------------------------------
    -- Recreate NOT deleted to inherit latest xref
    --------------------------------------------------------
    execute 'drop table if exists xref_pel_not_deleted';
    -- execute 'create table if not exists xref_pel_not_deleted (like xref including defaults including constraints including indexes including comments)';
    execute 'create table if not exists xref_pel_not_deleted (like xref including defaults including constraints including comments including generated)';

    --------------------------------------------------------
    --  Drop index if exists on NOT deleted table
    --------------------------------------------------------
    execute 'drop index if exists xref_pel_not_deleted$ac';
    execute 'drop index if exists xref_pel_not_deleted$created';
    execute 'drop index if exists xref_pel_not_deleted$upi';
    execute 'drop index if exists xref_pel_not_deleted$taxid';
    --------------------------------------------------------
    --  Drop FK constraints on NOT deleted table
    --------------------------------------------------------
    execute 'alter table xref_pel_not_deleted drop constraint if exists xref_pel_not_deleted_fk1';
    execute 'alter table xref_pel_not_deleted drop constraint if exists xref_pel_not_deleted_fk2';
    execute 'alter table xref_pel_not_deleted drop constraint if exists xref_pel_not_deleted_fk3';
    execute 'alter table xref_pel_not_deleted drop constraint if exists xref_pel_not_deleted_fk4';

    -- code here to raise notice/exception in case of any new index/constraint not managed by this procedure for both tables

  END;

$function$

