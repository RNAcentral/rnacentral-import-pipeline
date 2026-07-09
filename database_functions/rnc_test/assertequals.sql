CREATE OR REPLACE FUNCTION rnc_test.assertequals(proc text, expected bigint, actual bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    IF NOT expected = actual THEN
      RAISE NOTICE 'Warning! % expected %, got %', proc, expected, actual;
    ELSE
      RAISE NOTICE '% ok', proc;
    END IF;
  END;

$function$

