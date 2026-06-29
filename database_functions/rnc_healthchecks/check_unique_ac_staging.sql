CREATE OR REPLACE FUNCTION rnc_healthchecks.check_unique_ac_staging()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_ac_all bigint;
    v_ac_distinct bigint;

BEGIN
    SELECT count(*) into v_ac_all FROM load_rnacentral_all;
    SELECT count(DISTINCT ac) into v_ac_distinct FROM load_rnacentral_all;

    IF v_ac_all != v_ac_distinct THEN

      RAISE NOTICE 'not ok ... check_unique_ac_staging';
    ELSE
      RAISE NOTICE 'ok ... check_unique_ac_staging';
    end if;

  END;

$function$

