CREATE OR REPLACE FUNCTION rnc_test.import_staging_data(p_test_id bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    -- do not insert anything for this test.
    -- new release contains no records.
    IF p_test_id = 18 THEN
      RETURN;
    END IF;

    INSERT INTO
      RNACEN.load_rnacentral_all(
        crc64,
        len,
        seq_short,
        seq_long,
        ac,
        VERSION,
        taxid,
        md5,
        database
      )
      VALUES(current_setting('RNC_TEST.v_crc')::CrcList(p_test_id), -- crc64
             current_setting('RNC_TEST.SEQ_LENGTH')::bigint,       -- length
             current_setting('RNC_TEST.v_seq')::SeqList(p_test_id), -- seq
             NULL,             -- long seq
             current_setting('RNC_TEST.v_acc')::AccList(p_test_id), -- accession
             current_setting('RNC_TEST.v_ver')::VerList(p_test_id), -- version
             current_setting('RNC_TEST.v_tax')::TaxList(p_test_id), -- taxid
             current_setting('RNC_TEST.v_md5')::Md5List(p_test_id), -- md5
             current_setting('RNC_TEST.DEFAULT_DB_NAME')::varchar(3)   -- default database label
             );
    --COMMIT;
  END;

$function$

