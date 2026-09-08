CREATE OR REPLACE FUNCTION rnc_load_xref.do_pel_exchange(p_in_db_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE
    v_parent_partition text;
    v_deleted_partition text := 'xref_p' || p_in_db_id || '_deleted';
    v_not_deleted_partition text := 'xref_p' || p_in_db_id || '_not_deleted';
    v_constraint_name text;
BEGIN
    --------------------------------------------------------
    --  drop any existing old partition tables
    --------------------------------------------------------
    execute 'drop table if exists xref_p' || p_in_db_id || '_deleted_old';
    execute 'drop table if exists xref_p' || p_in_db_id || '_not_deleted_old';

    --------------------------------------------------------
    --  rename indexes on current partition tables as old
    --------------------------------------------------------
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$id rename to xref_p' || p_in_db_id || '_deleted$id_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$dbid rename to xref_p' || p_in_db_id || '_deleted$dbid_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$ac rename to xref_p' || p_in_db_id || '_deleted$ac_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$created rename to xref_p' || p_in_db_id || '_deleted$created_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$last rename to xref_p' || p_in_db_id || '_deleted$last_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$upi rename to xref_p' || p_in_db_id || '_deleted$upi_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_deleted$taxid rename to xref_p' || p_in_db_id || '_deleted$taxid_old';

    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$id rename to xref_p' || p_in_db_id || '_not_deleted$id_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$dbid rename to xref_p' || p_in_db_id || '_not_deleted$dbid_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$ac rename to xref_p' || p_in_db_id || '_not_deleted$ac_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$created rename to xref_p' || p_in_db_id || '_not_deleted$created_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$last rename to xref_p' || p_in_db_id || '_not_deleted$last_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$upi rename to xref_p' || p_in_db_id || '_not_deleted$upi_old';
    execute 'alter index if exists xref_p' || p_in_db_id || '_not_deleted$taxid rename to xref_p' || p_in_db_id || '_not_deleted$taxid_old';

    --------------------------------------------------------
    --  rename FK constraints on current partition tables as old
    --------------------------------------------------------
    FOREACH v_constraint_name IN ARRAY ARRAY[
        'xref_p' || p_in_db_id || '_deleted_fk1',
        'xref_p' || p_in_db_id || '_deleted_fk2',
        'xref_p' || p_in_db_id || '_deleted_fk3',
        'xref_p' || p_in_db_id || '_deleted_fk4'
    ]
    LOOP
        IF EXISTS (
            SELECT 1
            FROM pg_constraint c
            JOIN pg_class r ON r.oid = c.conrelid
            JOIN pg_namespace n ON n.oid = r.relnamespace
            WHERE n.nspname = current_schema()
              AND r.relname = v_deleted_partition
              AND c.conname = v_constraint_name
        ) THEN
            execute format(
                'alter table %I rename constraint %I to %I',
                v_deleted_partition,
                v_constraint_name,
                v_constraint_name || '_old'
            );
        END IF;
    END LOOP;

    FOREACH v_constraint_name IN ARRAY ARRAY[
        'xref_p' || p_in_db_id || '_not_deleted_fk1',
        'xref_p' || p_in_db_id || '_not_deleted_fk2',
        'xref_p' || p_in_db_id || '_not_deleted_fk3',
        'xref_p' || p_in_db_id || '_not_deleted_fk4'
    ]
    LOOP
        IF EXISTS (
            SELECT 1
            FROM pg_constraint c
            JOIN pg_class r ON r.oid = c.conrelid
            JOIN pg_namespace n ON n.oid = r.relnamespace
            WHERE n.nspname = current_schema()
              AND r.relname = v_not_deleted_partition
              AND c.conname = v_constraint_name
        ) THEN
            execute format(
                'alter table %I rename constraint %I to %I',
                v_not_deleted_partition,
                v_constraint_name,
                v_constraint_name || '_old'
            );
        END IF;
    END LOOP;

    --------------------------------------------------------
    -- find the actual partition parent for the current dbid branch
    --------------------------------------------------------
    select inhparent::regclass::text
      into v_parent_partition
      from pg_inherits
     where inhrelid = to_regclass(v_deleted_partition);

    if v_parent_partition is null then
        raise exception 'Could not determine xref partition parent for dbid %', p_in_db_id;
    end if;

    --------------------------------------------------------
    -- detach current partitions before renaming them away
    --------------------------------------------------------
    execute format('alter table %s detach partition %I', v_parent_partition, v_deleted_partition);
    execute format('alter table %s detach partition %I', v_parent_partition, v_not_deleted_partition);

    --------------------------------------------------------
    --  Rename current partition tables as old
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted rename to xref_p' || p_in_db_id || '_deleted_old';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted rename to xref_p' || p_in_db_id || '_not_deleted_old';

    --------------------------------------------------------
    --  Rename pel tables according to naming conventions for partition tables of xref
    --------------------------------------------------------
    execute 'alter table xref_pel_deleted rename to xref_p' || p_in_db_id || '_deleted';
    execute 'alter table xref_pel_not_deleted rename to xref_p' || p_in_db_id || '_not_deleted';

    --------------------------------------------------------
    --  Add check constraints to accept only data complying with partitions definition
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_check' ||
            ' check ((dbid = ' || p_in_db_id || ') AND (deleted = ''Y''))';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_check' ||
            ' check ((dbid = ' || p_in_db_id || ') AND (deleted = ''N''))';

    --------------------------------------------------------
    --  Creating indexes on deleted partition
    --------------------------------------------------------
    execute 'create unique index xref_p' || p_in_db_id || '_deleted$id on xref_p' || p_in_db_id || '_deleted (id)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$dbid on xref_p' || p_in_db_id || '_deleted (dbid)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$ac on xref_p' || p_in_db_id || '_deleted (ac)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$created on xref_p' || p_in_db_id || '_deleted (created)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$last on xref_p' || p_in_db_id || '_deleted (last)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$upi on xref_p' || p_in_db_id || '_deleted (urs)';
    execute 'create index xref_p' || p_in_db_id || '_deleted$taxid on xref_p' || p_in_db_id || '_deleted (taxid)';

    --------------------------------------------------------
    --  Creating FK constraints on deleted partition
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk1' ||
            ' foreign key(created) references rnc_release (id)';
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk2' ||
            ' foreign key(dbid) references rnc_database (id)';
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk3' ||
            ' foreign key(last) references rnc_release (id)';
    -- urs -> rna is validated by construction (every urs written here came from
    -- comparable_prot_upi/prot_upi, which already exist in rna). Added NOT VALID to skip
    -- the in-transaction full-partition scan + probes into the 35GB rna(urs) index, which
    -- dominates the per-database load. run() runs VALIDATE CONSTRAINT after the load loop
    -- (ShareUpdateExclusive, off the critical path) to mark it valid and catch any violation.
    execute 'alter table xref_p' || p_in_db_id || '_deleted add constraint xref_p' || p_in_db_id || '_deleted_fk4' ||
            ' foreign key(urs) references rna (urs) not valid';

    --------------------------------------------------------
    --  Creating indexes on NOT deleted partition
    --------------------------------------------------------
    execute 'create unique index xref_p' || p_in_db_id || '_not_deleted$id on xref_p' || p_in_db_id || '_not_deleted (id)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$dbid on xref_p' || p_in_db_id || '_not_deleted (dbid)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$ac on xref_p' || p_in_db_id || '_not_deleted (ac)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$created on xref_p' || p_in_db_id || '_not_deleted (created)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$last on xref_p' || p_in_db_id || '_not_deleted (last)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$upi on xref_p' || p_in_db_id || '_not_deleted (urs)';
    execute 'create index xref_p' || p_in_db_id || '_not_deleted$taxid on xref_p' || p_in_db_id || '_not_deleted (taxid)';

    --------------------------------------------------------
    --  Creating FK constraints on NOT deleted partition
    --------------------------------------------------------
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk1' ||
            ' foreign key(created) references rnc_release (id)';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk2' ||
            ' foreign key(dbid) references rnc_database (id)';
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk3' ||
            ' foreign key(last) references rnc_release (id)';
    -- See note on _deleted_fk4 above: NOT VALID here, validated off the critical path in run().
    execute 'alter table xref_p' || p_in_db_id || '_not_deleted add constraint xref_p' || p_in_db_id || '_not_deleted_fk4' ||
            ' foreign key(urs) references rna (urs) not valid';

    --------------------------------------------------------
    --  attach fresh partitions back to the correct parent
    --------------------------------------------------------
    execute format('alter table %s attach partition %I for values in (''Y'')', v_parent_partition, v_deleted_partition);
    execute format('alter table %s attach partition %I for values in (''N'')', v_parent_partition, v_not_deleted_partition);

    --------------------------------------------------------
    --   gathering statistics
    --------------------------------------------------------
    execute 'analyze xref_p' || p_in_db_id || '_deleted';
    execute 'analyze xref_p' || p_in_db_id || '_not_deleted';
END;
$function$
