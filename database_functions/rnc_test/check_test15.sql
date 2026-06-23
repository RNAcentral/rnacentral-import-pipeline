CREATE OR REPLACE FUNCTION rnc_test.check_test15()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    perform check_test9;

  END;

$function$

