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

    perform rnc_load_xref.load_upi_max_versions_table(p_in_dbid);
    perform rnc_load_xref.load_max_versions_table();

    perform rnc_load_xref_incremental.incremental_new_versions(p_previous_release);
    perform rnc_load_xref_incremental.incremental_new_accessions(p_in_dbid);

    perform rnc_load_xref_incremental.incremental_refresh(p_previous_release);
    perform rnc_load_xref_incremental.incremental_retire_changed(p_previous_release);
    perform rnc_load_xref_incremental.incremental_retire_explicit(p_in_dbid, p_previous_release);
END;
$function$
