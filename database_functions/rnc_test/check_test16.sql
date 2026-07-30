CREATE OR REPLACE FUNCTION rnc_test.check_test16()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN


    perform check_test12;

  END;

$function$

