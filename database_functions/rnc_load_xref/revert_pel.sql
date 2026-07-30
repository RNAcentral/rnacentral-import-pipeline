CREATE OR REPLACE FUNCTION rnc_load_xref.revert_pel(p_in_db_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    --------------------------------------------------------
    --  Swap inherit clause to disable/enable data selection in query submitted
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted_old inherit xref';    
    execute 'alter table if exists xref_p' || p_in_db_id || '_deleted no inherit xref';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted_old inherit xref';    
    execute 'alter table if exists xref_p' || p_in_db_id || '_not_deleted no inherit xref';

    --------------------------------------------------------
    --  drop partition tables not inheriting anymore
    --------------------------------------------------------
    execute 'drop table if exists xref_p' || p_in_db_id || '_deleted';
    execute 'drop table if exists xref_p' || p_in_db_id || '_not_deleted';
    
    --------------------------------------------------------
    --  Rename old partition tables back to current
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted_old rename to xref_p' || p_in_db_id || '_deleted';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted_old rename to xref_p' || p_in_db_id || '_not_deleted';

    --------------------------------------------------------
    --  rename indexes on current partition tables as old
    --------------------------------------------------------
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$id_old rename to xref_p' || p_in_db_id || '_deleted$id';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$ac_old rename to xref_p' || p_in_db_id || '_deleted$ac';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$created_old rename to xref_p' || p_in_db_id || '_deleted$created';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$upi_old rename to xref_p' || p_in_db_id || '_deleted$upi';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$taxid_old rename to xref_p' || p_in_db_id || '_deleted$taxid';

    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$id_old rename to xref_p' || p_in_db_id || '_not_deleted$id';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$ac_old rename to xref_p' || p_in_db_id || '_not_deleted$ac';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$created_old rename to xref_p' || p_in_db_id || '_not_deleted$created';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$upi_old rename to xref_p' || p_in_db_id || '_not_deleted$upi';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$taxid_old rename to xref_p' || p_in_db_id || '_not_deleted$taxid';

    --------------------------------------------------------
    --  rename FK constraints on current partition tables as old
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk1_old to xref_p' || p_in_db_id || '_deleted_fk1';
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk2_old to xref_p' || p_in_db_id || '_deleted_fk2';
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk3_old to xref_p' || p_in_db_id || '_deleted_fk3';
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename constraint xref_p' || p_in_db_id || '_deleted_fk4_old to xref_p' || p_in_db_id || '_deleted_fk4';

    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk1_old to xref_p' || p_in_db_id || '_not_deleted_fk1';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk2_old to xref_p' || p_in_db_id || '_not_deleted_fk2';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk3_old to xref_p' || p_in_db_id || '_not_deleted_fk3';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename constraint xref_p' || p_in_db_id || '_not_deleted_fk4_old to xref_p' || p_in_db_id || '_not_deleted_fk4';

END;

$function$

