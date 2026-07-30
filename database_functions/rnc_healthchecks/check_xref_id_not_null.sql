CREATE OR REPLACE FUNCTION rnc_healthchecks.check_xref_id_not_null()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_count bigint;

BEGIN
    SELECT count(*) INTO v_count FROM xref WHERE id IS NULL;
    perform rnc_healthchecks.report_results(v_count, 'check_xref_id_not_null');
  END;

$function$

