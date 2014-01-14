set define off

create or replace PACKAGE RNC_TEST AS

  /*
    Package for testing the RNAcentral data import pipeline.
  */

  PROCEDURE run_tests;

END RNC_TEST;
/
create or replace PACKAGE BODY RNC_TEST AS





  /*
    Warning: this package truncates all critical tables
    in order to run tests in a clean environment.
    The tests should never be run in production.

    Input records in the staging table contain 4 primary fields:
    sequence, accession, version, and taxid.
    Other fields are derived from these 4.

    Therefore, given an input record in one release, there are 16 ways
    in which this record can change in the next release.

    These tests are designed to cover all 16 possibilities.


    Find additional documentation at http://goo.gl/U1bWK

    Currently not covered:
      - sequences longer than 4000 nts
      - multiple preceding xrefs
      - incremental releases
  */

  /*
    Global constants

  */
  SEQ_LENGTH CONSTANT NUMBER := 5;
  FIRST_UPI  CONSTANT VARCHAR(13) := 'UPI0000000001';

  DEFAULT_DB_NAME CONSTANT VARCHAR(3) := 'ENA';

  /*
    Test parameters (available globally).

    The first entry in all arrays is the default record that everything else
    is compared to.
  */


  -- define array types
  TYPE SeqList IS VARRAY(17) OF VARCHAR(10); -- first entry + 16 test cases
  TYPE AccList IS VARRAY(17) OF VARCHAR(4);
  TYPE VerList IS VARRAY(17) OF NUMBER;

  TYPE TaxList IS VARRAY(17) OF NUMBER;
  TYPE CrcList IS VARRAY(17) OF VARCHAR(16);
  TYPE Md5List IS VARRAY(17) OF VARCHAR(32);

   -- first 9 sequences are identical, the next 8 are different, 17 total
  v_seq SeqList := SeqList('AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA',
                           'GGGGG','GGGGG','GGGGG','GGGGG','GGGGG','GGGGG','GGGGG','GGGGG');
  v_acc AccList := AccList('id1','id1','id2','id1','id1','id2','id2','id1','id2','id2','id2','id2','id1','id1','id1','id2','id1');

  v_ver VerList := VerList(   1 ,   1 ,   1 ,   2 ,   1 ,   2 ,   1 ,   2 ,   2 ,   2 ,   2 ,   1 ,   2 ,   1 ,   2 ,   1 ,   1);
  v_tax TaxList := TaxList(   1 ,   1 ,   1 ,   1 ,   2 ,   1 ,   2 ,   2 ,   2 ,   2 ,   1 ,   2 ,   2 ,   2 ,   1 ,   1 ,   1);
  v_crc CrcList := CrcList('0000000000000000','0000000000000000','0000000000000000',
                           '0000000000000000','0000000000000000','0000000000000000',
                           '0000000000000000','0000000000000000','0000000000000000',

                           '1111111111111111','1111111111111111','1111111111111111','1111111111111111',
                           '1111111111111111','1111111111111111','1111111111111111','1111111111111111'
                           );
  v_md5 Md5List := Md5List('00000000000000000000000000000000','00000000000000000000000000000000','00000000000000000000000000000000',
                           '00000000000000000000000000000000','00000000000000000000000000000000','00000000000000000000000000000000',
                           '00000000000000000000000000000000','00000000000000000000000000000000','00000000000000000000000000000000',
                           '11111111111111111111111111111111','11111111111111111111111111111111','11111111111111111111111111111111','11111111111111111111111111111111',

                           '11111111111111111111111111111111','11111111111111111111111111111111','11111111111111111111111111111111','11111111111111111111111111111111'
                           );

  /* Insert the first record into the RNA table. */
  PROCEDURE initialize_rna_table AS
  BEGIN

    INSERT
      INTO rnacen.rna
        (ID,
         UPI,
         TIMESTAMP,
         USERSTAMP,

         CRC64,
         LEN,
         SEQ_SHORT,
         SEQ_LONG,
         MD5)
      VALUES(1,                 -- id
             FIRST_UPI,         -- UPI

             CURRENT_TIMESTAMP, -- timestamp
             USER,              -- user stamp
             v_crc(1),          -- crc64
             SEQ_LENGTH,        -- length
             v_seq(1),          -- sequence

             NULL,              -- long sequence
             v_md5(1)           -- md5
             );
    COMMIT;
  END initialize_rna_table;

  /* Insert the first record into the XREF table. */
  PROCEDURE initialize_xref_table AS

  BEGIN
    INSERT
      INTO RNACEN.xref
      (

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

      VALUES(
             1,                 -- dbid
             1,                 -- first release
             1,                 -- last release
             FIRST_UPI,         -- UPI
             1,                 -- version_I
             v_ver(1),          -- version
             'N',               -- deleted
             CURRENT_TIMESTAMP, -- timestamp
             USER,              -- userstamp

             v_acc(1),          -- accession
             v_tax(1)           -- taxid

             );
    COMMIT;
  END initialize_xref_table;

  /*
    Insert two release records. The first is the release corresponding
    to the initial entry in xref and in rna tables. The second release
    is the new release for the test entry.
  */
  PROCEDURE initialize_releases AS
  BEGIN

    -- add first release, set release status to "Done"

    INSERT
      INTO RNACEN.rnc_release
      VALUES(1,                 -- id
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
      VALUES(2,1,CURRENT_TIMESTAMP,'F','L',CURRENT_TIMESTAMP,USER,NULL,'N');

    COMMIT;
  END initialize_releases;

  /*
    Set up initial values
  */
  PROCEDURE setup AS


  BEGIN

    initialize_releases;
    initialize_rna_table;
    initialize_xref_table;

  END setup;

  /*
    Import new data into the staging table load_rnacentral
  */
  PROCEDURE import_staging_data(

    p_test_id NUMBER

  ) AS
  BEGIN

    -- do not insert anything for this test.
    -- new release contains no records.
    IF p_test_id = 18 THEN
      RETURN;
    END IF;

    INSERT INTO
      RNACEN.load_rnacentral_all

      (
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
      VALUES(v_crc(p_test_id), -- crc64

             SEQ_LENGTH,       -- length
             v_seq(p_test_id), -- seq
             NULL,             -- long seq

             v_acc(p_test_id), -- accession
             v_ver(p_test_id), -- version
             v_tax(p_test_id), -- taxid
             v_md5(p_test_id), -- md5
             DEFAULT_DB_NAME   -- default database label
             );
    COMMIT;
  END import_staging_data;


  /*
    Truncate all test tables
  */
  PROCEDURE teardown AS

  BEGIN

    EXECUTE IMMEDIATE 'TRUNCATE TABLE load_rnacentral_all DROP STORAGE';

    EXECUTE IMMEDIATE 'ALTER TABLE xref DISABLE CONSTRAINT FK_XREF$CREATED';
    EXECUTE IMMEDIATE 'ALTER TABLE xref DISABLE CONSTRAINT FK_XREF$DBID';
    EXECUTE IMMEDIATE 'ALTER TABLE xref DISABLE CONSTRAINT FK_XREF$LAST';
    EXECUTE IMMEDIATE 'ALTER TABLE xref DISABLE CONSTRAINT FK_XREF$UPI';


    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED DISABLE CONSTRAINT XREF_PEL_DELETED_FK1';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED DISABLE CONSTRAINT XREF_PEL_DELETED_FK2';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED DISABLE CONSTRAINT XREF_PEL_DELETED_FK3';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED DISABLE CONSTRAINT XREF_PEL_DELETED_FK4';

    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED DISABLE CONSTRAINT XREF_PEL_NOT_DELETED_FK1';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED DISABLE CONSTRAINT XREF_PEL_NOT_DELETED_FK2';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED DISABLE CONSTRAINT XREF_PEL_NOT_DELETED_FK3';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED DISABLE CONSTRAINT XREF_PEL_NOT_DELETED_FK4';

    EXECUTE IMMEDIATE 'ALTER TABLE RNACEN.rnc_release DISABLE CONSTRAINT RNC_RELEASE$ID$PK';


    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rnc_release DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rna DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.xref DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.load_rnacentral DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_DELETED DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_NOT_DELETED DROP STORAGE';


    EXECUTE IMMEDIATE 'ALTER TABLE RNACEN.rnc_release ENABLE CONSTRAINT RNC_RELEASE$ID$PK';

    EXECUTE IMMEDIATE 'ALTER TABLE xref ENABLE CONSTRAINT FK_XREF$CREATED';
    EXECUTE IMMEDIATE 'ALTER TABLE xref ENABLE CONSTRAINT FK_XREF$DBID';
    EXECUTE IMMEDIATE 'ALTER TABLE xref ENABLE CONSTRAINT FK_XREF$LAST';

    EXECUTE IMMEDIATE 'ALTER TABLE xref ENABLE CONSTRAINT FK_XREF$UPI';

    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED ENABLE CONSTRAINT XREF_PEL_DELETED_FK1';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED ENABLE CONSTRAINT XREF_PEL_DELETED_FK2';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED ENABLE CONSTRAINT XREF_PEL_DELETED_FK3';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_DELETED ENABLE CONSTRAINT XREF_PEL_DELETED_FK4';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED ENABLE novalidate CONSTRAINT XREF_PEL_NOT_DELETED_FK1';

    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED ENABLE novalidate CONSTRAINT XREF_PEL_NOT_DELETED_FK2';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED ENABLE novalidate CONSTRAINT XREF_PEL_NOT_DELETED_FK3';
    EXECUTE IMMEDIATE 'ALTER TABLE XREF_PEL_NOT_DELETED ENABLE novalidate CONSTRAINT XREF_PEL_NOT_DELETED_FK4';

  END teardown;


  /*
    Utility function for unit testing.
    Example: assertEquals(2, 2, 'function_name');
    Adopted from https://github.com/uzilan/asserts-package
  */
  PROCEDURE assertEquals(
    proc     IN VARCHAR2,

    expected IN NUMBER,
    actual   IN NUMBER)
  IS
  BEGIN

    IF NOT expected = actual THEN
      DBMS_OUTPUT.put_line('Warning! ' || proc || ' expected ' || expected || ', got ' || actual);
    ELSE
      DBMS_OUTPUT.put_line(proc || ' ok');
    END IF;
  end assertEquals;

  /*
  *********

  * TESTS *
  *********
  */


  /*
    Input: identical entry.
    Expected action: increment Last seen field in the existing xref.
  */
  PROCEDURE check_test1 AS
    l_count NUMBER;
    l_test_id VARCHAR(20);
  BEGIN


    l_test_id := 'check_test1';


    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 1 AND last = 2;
    assertEquals(l_test_id, 1, l_count);

  END check_test1;



  /*
    Input: different accession number.
    Expected action: retire existing xref, create a new xref with new accession.
  */
  PROCEDURE check_test2 AS
    l_count NUMBER;
    l_test_id VARCHAR(20);
  BEGIN

    l_test_id := 'check_test2';

    SELECT count(*) INTO l_count FROM rna;


    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1 AND ac = 'id1';
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 2 AND LAST = 2 AND ac = 'id2';
    assertEquals(l_test_id, 1, l_count);



  END check_test2;

  /*
    Input: different version number.
    Expected action: retire existing xref, create a new xref with updated version.
  */
  PROCEDURE check_test3 AS
    l_count NUMBER;
    l_test_id VARCHAR(20);
  BEGIN

    l_test_id := 'check_test3';



    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1 AND version = 1;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref

    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 2 AND LAST = 2 AND version = 2;
    assertEquals(l_test_id, 1, l_count);


  END check_test3;

  /*
    Input: different taxid.
    Expected action: retire existing xref, create a new xref with updated taxid.
    Real action: overriding taxid in the existing xref, no change to version_I
    This test can fail.
  */
  PROCEDURE check_test4 AS

    l_count NUMBER;
    l_test_id VARCHAR(20);
  BEGIN


    l_test_id := 'check_test4';

    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);


    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1 AND taxid = 1;
    assertEquals(l_test_id, 1, l_count);


    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 2 AND LAST = 2 AND taxid = 2;
    assertEquals(l_test_id, 1, l_count);

  END check_test4;

  /*
    Input: Same seq and taxid, but diff accession and version.

    Expected action: retire existing xref, create a new xref with updated version and accession.
  */
  PROCEDURE check_test5 AS
    l_count NUMBER;
    l_test_id VARCHAR(20);

  BEGIN

    l_test_id := 'check_test5';

    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);


    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1 AND ac = 'id1' AND version = 1;
    assertEquals(l_test_id, 1, l_count);


    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 2 AND LAST = 2 AND ac = 'id2' AND version = 2;
    assertEquals(l_test_id, 1, l_count);

  END check_test5;


  /*
    Input: Same seq and version, but new acc and taxid
    Expected action: retire existing xref, create a new one with updated acc and taxid
  */
  PROCEDURE check_test6 AS
    l_count NUMBER;

    l_test_id VARCHAR(20);
  BEGIN

    l_test_id := 'check_test6';


    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1 AND ac = 'id1' AND taxid = 1;

    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 2 AND LAST = 2 AND ac = 'id2' AND taxid = 2;

    assertEquals(l_test_id, 1, l_count);

  END check_test6;

  /*
    Input: Same seq and acc, but new version and taxid
    Expected action: retire existing xref, create a new one with updated version and taxid
  */
  PROCEDURE check_test7 AS

    l_count NUMBER;
    l_test_id VARCHAR(20);
  BEGIN


    l_test_id := 'check_test7';

    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref

    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1 AND version = 1 AND taxid = 1;
    assertEquals(l_test_id, 1, l_count);


    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 2 AND LAST = 2 AND version = 2 AND taxid = 2;
    assertEquals(l_test_id, 1, l_count);

  END check_test7;

  /*
    Input: Same seq, but diff accession, ver, taxid
    Expected action: retire existing xref, create a new one with updated acc, version, and taxid
  */

  PROCEDURE check_test8 AS

    l_count NUMBER;
    l_test_id VARCHAR(20);
  BEGIN

    l_test_id := 'check_test8';

    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);



    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1 AND version = 1 AND taxid = 1 AND ac = 'id1';
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 2 AND LAST = 2 AND version = 2 AND taxid = 2 AND ac = 'id2';
    assertEquals(l_test_id, 1, l_count);

  END check_test8;

  /*
    Tests 8-16 are for entries with different sequences.
  /*



  /*
    Input: new seq, ac, ver, taxid
    Expected action: retire existing xref because it's no longer present
    in the release, add a new sequence, add a new xref.

    This test is also run on test entries 10, 11 and 15.
  */
  PROCEDURE check_test9 AS
    l_count NUMBER;
    l_test_id VARCHAR(30);
  BEGIN



    l_test_id := 'check_test_9,10,11,15';

    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1;
    assertEquals(l_test_id, 1, l_count);


    SELECT count(*) INTO l_count FROM xref

    WHERE upi != FIRST_UPI AND deleted = 'N' AND created = 2;
    assertEquals(l_test_id, 1, l_count);

  END check_test9;

  /*
    Input: Only tax id is the same.
    Expected action: New seq and xref, retire existing xref.
  */
  PROCEDURE check_test10 AS

  BEGIN

    check_test9;


  END check_test10;

  /*
    Input: Only version is the same.
    Expected action: New seq and xref, retire existing xref.
  */
  PROCEDURE check_test11 AS
  BEGIN


    check_test9;

  END check_test11;


  /*
    Input: Same acc, but new seq, ver, taxid.
    Expected action: New seq and xref, retire existing xref. Increment VERSION_I.
    Applies to tests 12,13,14,16
  */
  PROCEDURE check_test12 AS
    l_count NUMBER;

    l_test_id VARCHAR(30);
  BEGIN

    l_test_id := 'check_test12_12,13,14,16';


    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 2, l_count);

    SELECT count(*) INTO l_count FROM xref

    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND version_I = 1;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi != FIRST_UPI AND deleted = 'N' AND created = 2 AND version_I = 2;
    assertEquals(l_test_id, 1, l_count);


  END check_test12;

  /*
    Input: Same acc and ver, but new seq and taxid
    Expected action: New seq and xref, retire existing xref. Increment VERSION_I.

  */
  PROCEDURE check_test13 AS
  BEGIN

    check_test12;

  END check_test13;


  /*
    Input: Same acc and taxid, but updated seq and version
    Expected action: New seq and xref, retire existing xref. Increment VERSION_I.
  */

  PROCEDURE check_test14 AS
  BEGIN

    check_test12;

  END check_test14;

  /*

    Input: Only version and tax id are the same
    Expected action: New seq and xref, retire existing xref.
  */
  PROCEDURE check_test15 AS

  BEGIN

    check_test9;

  END check_test15;

  /*
    Input: same ac, ver, taxid, but updated seq
    Expected action: New seq and xref, retire existing xref. Increment VERSION_I.

  */
  PROCEDURE check_test16 AS
  BEGIN


    check_test12;

  END check_test16;

  /*
    Input: new release contains no entries.
    Expected action: retire existing xref.
  */
  PROCEDURE check_test17 AS

    l_count NUMBER;
    l_test_id VARCHAR(20);

  BEGIN

    l_test_id := 'check_test17';

    SELECT count(*) INTO l_count FROM rna;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals(l_test_id, 1, l_count);

    SELECT count(*) INTO l_count FROM xref

    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1;

    assertEquals(l_test_id, 1, l_count);

  END check_test17;

  /*
    Check results
  */
  PROCEDURE check_result(
    p_test_id IN NUMBER
  ) AS
  BEGIN



    -- dynamic function call turned out to be tricky to set up
    CASE p_test_id
      WHEN 1  THEN check_test1;
      WHEN 2  THEN check_test2;
      WHEN 3  THEN check_test3;
      WHEN 4  THEN check_test4;
      WHEN 5  THEN check_test5;
      WHEN 6  THEN check_test6;
      WHEN 7  THEN check_test7;
      WHEN 8  THEN check_test8;
      WHEN 9  THEN check_test9;
      WHEN 10 THEN check_test10;
      WHEN 11 THEN check_test11;


      WHEN 12 THEN check_test12;
      WHEN 13 THEN check_test13;
      WHEN 14 THEN check_test14;
      WHEN 15 THEN check_test15;
      WHEN 16 THEN check_test16;
      WHEN 17 THEN check_test17;
      ELSE
        DBMS_OUTPUT.put_line('Unknown test id ' || p_test_id);
    END CASE;

  END check_result;


  /*

    Run a single test specified in the parameter
  */
  PROCEDURE run_test(
    p_test_id NUMBER
  ) AS
  BEGIN

    -- initialize rna, xref, and release tables
    setup;

    -- insert new data into the staging table

    -- +1 accounts for the first (default) entry
    import_staging_data( p_test_id + 1 );


    -- run update
    RNC_UPDATE.new_update;

    -- check results
    check_result(p_test_id);

    -- clean up
    teardown;


  END run_test;

  /*

    Launch the test suite
  */
  PROCEDURE run_tests AS
    l_cntr NUMBER;
  BEGIN

    -- truncate all tables in the beginning
    teardown;


    -- 17 test cases
    FOR l_cntr IN 1..17
    LOOP
      DBMS_OUTPUT.put_line('Running test ' || l_cntr);

      run_test(l_cntr);
    END LOOP;

  END run_tests;


  /*
    Immediately truncate all tables

  */
  PROCEDURE truncate_all_tables AS
  begin
    teardown;
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.release_stats DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rnc_ac_info DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rnc_accessions DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rnc_composite_ids DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rnc_assembly DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rnc_references DROP STORAGE';
    DBMS_OUTPUT.put_line('All tables truncated');
  end;



END RNC_TEST;
/
set define on
