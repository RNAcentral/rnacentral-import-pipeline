CREATE OR REPLACE FUNCTION rnc_load_rna.store_new_sequences()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    INSERT INTO rna(
        id,
        upi,
        crc64,
        LEN,
        seq_short,
        seq_long,
        md5,
        TIMESTAMP,
        USERSTAMP
      )
    SELECT 
          prot_id,
          prot_upi,
          in_crc64,
          in_len,
          in_seq_short,
          in_seq_long,
          in_md5,
          clock_timestamp(),
          USER
    from (
        SELECT DISTINCT
          n.PROT_ID,
          n.PROT_UPI,
          l.in_crc64,
          l.in_len,
          l.in_seq_short,
          l.IN_SEQ_LONG,
          l.in_md5
        FROM
          load_retro_tmp l,
          load_md5_new_sequences n
        WHERE
          n.in_md5 = l.in_md5) alias1;

    --COMMIT;

  END;

  /*
    The main procedure for importing data from the staging table
    into the main RNA table. Responsible for assigning UPIs.
  */


$function$

