CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.load_xref_delta(p_previous_release bigint, p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    -- DELTA load: the input holds only new/changed records, so deletions come from an
    -- explicit list (load_deletions), not absence. Reuses four incremental steps and
    -- swaps retire_dropped for retire_explicit; absent rows stay active.
    -- See docs/incremental-parsing.md.

    RAISE NOTICE 'load_xref_delta 1/7: load_upi_max_versions_table';
    perform rnc_load_xref.load_upi_max_versions_table(p_in_dbid);
    RAISE NOTICE 'load_xref_delta 2/7: load_max_versions_table';
    perform rnc_load_xref.load_max_versions_table();

    RAISE NOTICE 'load_xref_delta 3/7: incremental_new_versions';
    perform rnc_load_xref_incremental.incremental_new_versions(p_in_dbid);
    RAISE NOTICE 'load_xref_delta 4/7: incremental_new_accessions';
    perform rnc_load_xref_incremental.incremental_new_accessions(p_in_dbid);

    RAISE NOTICE 'load_xref_delta 5/7: incremental_refresh';
    perform rnc_load_xref_incremental.incremental_refresh(p_in_dbid);
    RAISE NOTICE 'load_xref_delta 6/7: incremental_retire_changed';
    perform rnc_load_xref_incremental.incremental_retire_changed(p_in_dbid, p_previous_release);
    RAISE NOTICE 'load_xref_delta 7/7: incremental_retire_explicit';
    perform rnc_load_xref_incremental.incremental_retire_explicit(p_in_dbid, p_previous_release);
END;
$function$
