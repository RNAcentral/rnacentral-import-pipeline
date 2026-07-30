CREATE OR REPLACE FUNCTION rnc_test.initialize_xref_table()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN
    INSERT
      INTO RNACEN.xref(

        DBID,
        CREATED,
        LAST,
        UPI,
        VERSION_I,
        VERSION,
        DELETED,
        TIMESTAMP,
        USERSTAMP,

        AC,
        TAXID
      )

      VALUES (
             1,                 -- dbid
             1,                 -- first release
             1,                 -- last release
             current_setting('RNC_TEST.FIRST_UPI')::varchar(13),         -- UPI
             1,                 -- version_I
             current_setting('RNC_TEST.v_ver')::VerList(1),          -- version
             'N',               -- deleted
             CURRENT_TIMESTAMP, -- timestamp
             USER,              -- userstamp

             current_setting('RNC_TEST.v_acc')::AccList(1),          -- accession
             current_setting('RNC_TEST.v_tax')::TaxList(1)           -- taxid

             );
    --COMMIT;
  END;

$function$

