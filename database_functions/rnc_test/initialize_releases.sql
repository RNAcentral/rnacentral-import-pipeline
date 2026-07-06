CREATE OR REPLACE FUNCTION rnc_test.initialize_releases()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    -- add first release, set release status to "Done"

    INSERT
      INTO RNACEN.rnc_release
      VALUES (1,                 -- id
             1,                 -- dbid
             CURRENT_TIMESTAMP, -- release date
             'F',               -- release type
             'D',               -- release status
             CURRENT_TIMESTAMP, -- timestamp
             USER,              -- userstamp
             'test release 1',  -- description
             'N'                -- force load
             );



    -- add second release
    INSERT INTO
      RNACEN.rnc_release
      VALUES (2,1,CURRENT_TIMESTAMP,'F','L',CURRENT_TIMESTAMP,USER,NULL,'N');

    --COMMIT;
  END;

$function$

