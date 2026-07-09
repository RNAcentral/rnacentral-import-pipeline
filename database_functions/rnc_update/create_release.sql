CREATE OR REPLACE FUNCTION rnc_update.create_release(p_in_dbid bigint, p_release_type character)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    v_next_release numeric;

BEGIN

    SELECT count(*) + 1 INTO v_next_release FROM rnc_release;


    RAISE NOTICE 'Creating new release for database %', p_in_dbid;

    INSERT INTO rnc_release(ID,
       dbid,
       release_date,
       release_type,
       status,
       TIMESTAMP,
       userstamp,
       descr,
       force_load)
    VALUES
      (
      v_next_release,
      p_in_dbid ,
      date_trunc('day', clock_timestamp()),
      p_release_type,
      'L',
      clock_timestamp(),
      'auto',
      '',
      'N'
    );

  -- --COMMIT;

  END;

$function$

