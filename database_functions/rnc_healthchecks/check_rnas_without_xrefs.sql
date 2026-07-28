CREATE OR REPLACE FUNCTION rnc_healthchecks.check_rnas_without_xrefs()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_count bigint;

BEGIN
    SELECT count(*) INTO v_count FROM xref x
RIGHT OUTER JOIN rna r ON (x.urs = r.urs)
WHERE x.urs IS NULL;

    IF v_count != 0 THEN
      RAISE NOTICE 'not ok ... check_rnas_without_xrefs';
    else
      RAISE NOTICE 'ok ... check_rnas_without_xrefs';
    END IF;


  END;

$function$

