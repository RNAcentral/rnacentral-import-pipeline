CREATE OR REPLACE FUNCTION rnc_load_rna.load_md5_new_sequences_table()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    EXECUTE 'TRUNCATE TABLE load_md5_new_sequences';

    INSERT INTO RNACEN.load_md5_new_sequences
      (
        IN_MD5,
        PROT_ID,
        PROT_UPI
      )
    SELECT
      IN_MD5,
      nextval('seq_upi') PROT_ID,
      upi.getUpi(currval('seq_upi')) PROT_UPI
    FROM
      LOAD_MD5_STATS
    WHERE
      NOT (
        CNT_DST_SEQ_SHORT > 1 OR CNT_DST_SEQ_LONG > 1
      );

    --COMMIT;

    execute 'analyze load_md5_new_sequences';

  END;

$function$

