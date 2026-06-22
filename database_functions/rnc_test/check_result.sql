CREATE OR REPLACE FUNCTION rnc_test.check_result(p_test_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN



    -- dynamic function call turned out to be tricky to set up
    CASE p_test_id
      WHEN 1  THEN perform check_test1;
      WHEN 2  THEN perform check_test2;
      WHEN 3  THEN perform check_test3;
      WHEN 4  THEN perform check_test4;
      WHEN 5  THEN perform check_test5;
      WHEN 6  THEN perform check_test6;
      WHEN 7  THEN perform check_test7;
      WHEN 8  THEN perform check_test8;
      WHEN 9  THEN perform check_test9;
      WHEN 10 THEN perform check_test10;
      WHEN 11 THEN perform check_test11;
      WHEN 12 THEN perform check_test12;
      WHEN 13 THEN perform check_test13;
      WHEN 14 THEN perform check_test14;
      WHEN 15 THEN perform check_test15;
      WHEN 16 THEN perform check_test16;
      WHEN 17 THEN perform check_test17;
      ELSE
        RAISE NOTICE 'Unknown test id %', p_test_id;
    END CASE;

  END;

$function$

