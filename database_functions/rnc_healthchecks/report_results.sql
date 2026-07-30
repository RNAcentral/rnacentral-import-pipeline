CREATE OR REPLACE FUNCTION rnc_healthchecks.report_results(p_in_count bigint, p_in_msg text)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    IF p_in_count > 0 THEN
      RAISE NOTICE 'not ok ... %', p_in_msg;
    ELSE
      RAISE NOTICE 'ok ... %', p_in_msg;
    END IF;

  END;

$function$

