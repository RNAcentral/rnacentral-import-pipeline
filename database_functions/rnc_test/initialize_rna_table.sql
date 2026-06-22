CREATE OR REPLACE FUNCTION rnc_test.initialize_rna_table()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    INSERT
      INTO rnacen.rna(ID,
         UPI,
         TIMESTAMP,
         USERSTAMP,

         CRC64,
         LEN,
         SEQ_SHORT,
         SEQ_LONG,
         MD5)
      VALUES(1,                 -- id
             current_setting('RNC_TEST.FIRST_UPI')::varchar(13),         -- UPI

             CURRENT_TIMESTAMP, -- timestamp
             USER,              -- user stamp
             current_setting('RNC_TEST.v_crc')::CrcList(1),          -- crc64
             current_setting('RNC_TEST.SEQ_LENGTH')::bigint,        -- length
             current_setting('RNC_TEST.v_seq')::SeqList(1),          -- sequence

             NULL,              -- long sequence
             current_setting('RNC_TEST.v_md5')::Md5List(1)           -- md5
             );
    --COMMIT;
  END;

$function$

