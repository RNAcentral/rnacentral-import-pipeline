CREATE OR REPLACE FUNCTION rnc_update.load_release(p_in_dbid bigint, p_in_load_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_previous_release RNACEN.rnc_release.ID%TYPE;
    v_release_type RNACEN.rnc_release.release_type%TYPE;

BEGIN

    RAISE NOTICE 'Loading release: %', p_in_load_release;

    -- initial logging
    perform rnc_logging.log_release_start(p_in_dbid, p_in_load_release);

    -- load sequences
    perform rnc_load_rna.load_rna(p_in_dbid, p_in_load_release);

    -- load xrefs
    v_release_type := release.get_release_type(p_in_load_release);
    v_previous_release := release.get_previous_release(p_in_dbid, p_in_load_release);
    IF v_release_type = 'F' THEN
      perform rnc_load_xref.load_xref(v_previous_release, p_in_dbid);
    elsif v_release_type = 'I' THEN
      perform rnc_load_xref_incremental.load_xref_incremental(v_previous_release, p_in_dbid);
    END IF;

    -- mark release as done
    perform rnc_update.mark_as_done(p_in_dbid, p_in_load_release);

    -- final logging
    perform rnc_logging.log_release_end(p_in_dbid, p_in_load_release, v_previous_release);

    --COMMIT;

  END;

$function$

