CREATE OR REPLACE FUNCTION rnc_healthchecks.check_xrefs_without_lit_refs()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_count bigint;

BEGIN
    select count(distinct urs) INTO v_count
    FROM rnc_reference_map t2
RIGHT OUTER JOIN xref t1 ON (t2.ACCESSION = t1.ac)
WHERE t2.accession is null and deleted = 'N';

    -- some entries don't have literature refs in ENA
    IF v_count > 100 THEN
      RAISE NOTICE 'not ok ... check_xrefs_without_literature_refs';
    else
      RAISE NOTICE 'ok ... check_xrefs_without_literature_refs';
    END IF;

  END;

$function$

