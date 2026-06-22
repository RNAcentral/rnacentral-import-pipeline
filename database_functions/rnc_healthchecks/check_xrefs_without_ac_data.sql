CREATE OR REPLACE FUNCTION rnc_healthchecks.check_xrefs_without_ac_data()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_count bigint;

BEGIN
    SELECT count(*) INTO v_count
    FROM rnc_accessions t2
RIGHT OUTER JOIN xref t1 ON (t2.ACCESSION = t1.ac)
WHERE t2.accession is null;

    IF v_count > 0 THEN
      RAISE NOTICE 'not ok ... check_xrefs_without_ac_data';
    else
      RAISE NOTICE 'ok ... check_xrefs_without_ac_data';
    END IF;

  END;

$function$

