CREATE OR REPLACE FUNCTION rnc_load_rna.store_new_sequences()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    INSERT INTO rna(
        id,
        urs,
        crc64,
        LEN,
        seq_short,
        seq_long,
        md5,
        TIMESTAMP,
        USERSTAMP
      )
    SELECT DISTINCT ON (n.in_md5)
          n.prot_id,
          n.prot_upi,
          l.in_crc64,
          l.in_len,
          l.in_seq_short,
          l.in_seq_long,
          n.in_md5,
          clock_timestamp(),
          USER
    FROM load_md5_new_sequences n
    JOIN load_retro_tmp l ON l.in_md5 = n.in_md5;

  END;
$function$
