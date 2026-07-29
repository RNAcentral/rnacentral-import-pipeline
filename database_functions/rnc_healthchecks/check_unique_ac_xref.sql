CREATE OR REPLACE FUNCTION rnc_healthchecks.check_unique_ac_xref()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_ac_all bigint;

    v_ac_distinct bigint;

BEGIN
    SELECT count(*) INTO v_ac_all FROM xref WHERE deleted='N';
    SELECT count(DISTINCT ac) into v_ac_distinct FROM xref where deleted='N';

    IF v_ac_all != v_ac_distinct THEN
      RAISE NOTICE 'not ok ... check_unique_ac_xref';
    ELSE
      RAISE NOTICE 'ok ... check_unique_ac_xref';
    end if;

  END;

$function$

