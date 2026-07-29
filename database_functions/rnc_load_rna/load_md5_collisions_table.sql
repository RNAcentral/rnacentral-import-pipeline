CREATE OR REPLACE FUNCTION rnc_load_rna.load_md5_collisions_table()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    EXECUTE 'TRUNCATE TABLE load_md5_collisions';

    INSERT
    INTO
      RNACEN.LOAD_MD5_COLLISIONS
      (
        IN_MD5,
        CNT,
        CNT_DST_SEQ_SHORT,
        CNT_DST_SEQ_LONG
      )
    SELECT
      IN_MD5,
      cnt,
      CNT_DST_SEQ_SHORT,
      CNT_DST_SEQ_LONG
    FROM
      LOAD_MD5_STATS
    WHERE (
        CNT_DST_SEQ_SHORT > 1
      OR CNT_DST_SEQ_LONG > 1
      );

    --COMMIT;

  END;

$function$

