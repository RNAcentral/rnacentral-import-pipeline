CREATE OR REPLACE FUNCTION rnc_healthchecks.check_negative_taxid()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_xref bigint;
    v_staging bigint;


BEGIN
    SELECT count(*) INTO v_xref FROM xref WHERE taxid < 0;
    SELECT count(*) into v_staging FROM load_rnacentral_all WHERE taxid < 0;

    IF v_xref != v_staging OR v_xref != 0 THEN
      RAISE NOTICE 'not ok ... check_negative_taxid';
    ELSE
      RAISE NOTICE 'ok ... check_negative_taxid';
    end if;

  END;

$function$

