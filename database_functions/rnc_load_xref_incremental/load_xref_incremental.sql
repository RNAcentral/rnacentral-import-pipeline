CREATE OR REPLACE FUNCTION rnc_load_xref_incremental.load_xref_incremental(p_previous_release bigint, p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    /*
    -- this is an INCREMENTAL release
    -- this is NOT a FULL release
    -- better to AVOID Partition Exchange Loading

    -- updates last, deleted and taxid
    -- of existing xrefs
    */

    perform rnc_load_xref_incremental.incremental1(p_previous_release);

    perform rnc_load_xref_incremental.incremental2(p_previous_release);

    perform rnc_load_xref_incremental.incremental3();

  END;

$function$

