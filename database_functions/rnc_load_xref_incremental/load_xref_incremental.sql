CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.load_xref_incremental(p_previous_release bigint, p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    -- Incremental xref load: touch only new/changed rows, leave everything else (incl.
    -- deleted history) in place -- avoids the full PEL rewrite. Each step mirrors a
    -- populate_pel_tables* step, so the active set matches a full release. Order matters
    -- (snapshot, insert, then update in place); see docs/incremental-xref-loading.md.

    -- 1. snapshot (identical to the tables the full path builds)
    RAISE NOTICE 'load_xref_incremental 1/7: load_upi_max_versions_table';
    perform rnc_load_xref.load_upi_max_versions_table(p_in_dbid);
    RAISE NOTICE 'load_xref_incremental 2/7: load_max_versions_table';
    perform rnc_load_xref.load_max_versions_table();

    -- 2. inserts (case C, then cases D + Gap A)
    RAISE NOTICE 'load_xref_incremental 3/7: incremental_new_versions';
    perform rnc_load_xref_incremental.incremental_new_versions(p_in_dbid);
    RAISE NOTICE 'load_xref_incremental 4/7: incremental_new_accessions';
    perform rnc_load_xref_incremental.incremental_new_accessions(p_in_dbid);

    -- 3. in-place updates (case 1 refresh, case 2 retire-changed, case 5 retire-dropped)
    RAISE NOTICE 'load_xref_incremental 5/7: incremental_refresh';
    perform rnc_load_xref_incremental.incremental_refresh(p_in_dbid);
    RAISE NOTICE 'load_xref_incremental 6/7: incremental_retire_changed';
    perform rnc_load_xref_incremental.incremental_retire_changed(p_in_dbid, p_previous_release);
    RAISE NOTICE 'load_xref_incremental 7/7: incremental_retire_dropped';
    perform rnc_load_xref_incremental.incremental_retire_dropped(p_in_dbid, p_previous_release);
END;
$function$
