CREATE OR REPLACE FUNCTION rnc_update.move_staging_data(p_in_dbid bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    RAISE NOTICE 'move_staging_data: copying dbid % into load_rnacentral', p_in_dbid;

    EXECUTE 'TRUNCATE TABLE load_rnacentral';

    INSERT INTO load_rnacentral (
      SELECT DISTINCT CRC64,
             LEN,
             SEQ_SHORT,
             SEQ_LONG,
             DATABASE,
             AC,
             OPTIONAL_ID,
             VERSION,
             TAXID,
             MD5
      FROM load_rnacentral_all d1, rnc_database d2
      WHERE d1.DATABASE = d2.descr AND d2.ID = p_in_dbid
    );

    --COMMIT;

  END;

$function$
