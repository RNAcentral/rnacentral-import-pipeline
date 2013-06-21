create or replace
PACKAGE BODY RNC_TEST AS

  /*
    Warning: this package truncates all critical tables
    so that tests run in a clean environment. Never run in production.

    Find additional documentation at http://goo.gl/U1bWK.
  */


  /*
    Global constants
  */
  SEQ_LENGTH CONSTANT NUMBER := 5;
  FIRST_UPI  CONSTANT VARCHAR(13) := 'UPI0000000001';

  /*
    Test parameters (available globally).

    The first entry in all arrays is the default record that everything else
    is compared to.
  */
  TYPE SeqList IS VARRAY(17) OF VARCHAR(10); -- first entry + 16 test cases
  TYPE AccList IS VARRAY(17) OF VARCHAR(4);
  TYPE VerList IS VARRAY(17) OF NUMBER;
  TYPE TaxList IS VARRAY(17) OF NUMBER;
  TYPE CrcList IS VARRAY(17) OF VARCHAR(16);
  TYPE Md5List IS VARRAY(17) OF VARCHAR(32);

   -- 9 identical, then 8 different, 17 total
  v_seq SeqList := SeqList('AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA','AAAAA',
                           'GGGGG','GGGGG','GGGGG','GGGGG','GGGGG','GGGGG','GGGGG','GGGGG');
  v_acc AccList := AccList('id1','id1','id2','id1','id1','id2','id2','id1','id2','id2','id2','id2','id1','id1','id1','id2','id1');
  v_ver VerList := VerList(1,1,1,2,1,2,1,2,2,2,2,1,2,1,2,1,1);
  v_tax TaxList := TaxList(100,100,100,100,200,100,200,200,200,200,100,200,200,200,100,100,100);
  v_crc CrcList := CrcList('0000000000000000','0000000000000000','0000000000000000','0000000000000000',
                           '0000000000000000','0000000000000000','0000000000000000','0000000000000000',
                           '0000000000000001','0000000000000001','0000000000000001','0000000000000001',
                           '0000000000000001','0000000000000001','0000000000000001','0000000000000001',
                           '0000000000000001');
  v_md5 Md5List := Md5List('00000000000000000000000000000000','00000000000000000000000000000000','00000000000000000000000000000000','00000000000000000000000000000000',
                           '00000000000000000000000000000000','00000000000000000000000000000000','00000000000000000000000000000000','00000000000000000000000000000000',
                           '00000000000000000000000000000001','00000000000000000000000000000001','00000000000000000000000000000001','00000000000000000000000000000001',
                           '00000000000000000000000000000001','00000000000000000000000000000001','00000000000000000000000000000001','00000000000000000000000000000001',
                           '00000000000000000000000000000001');

  /* Insert the first record. */
  PROCEDURE initialize_rna_table AS
  BEGIN
    INSERT
      INTO rnacen.rna
      VALUES(1, -- id
             FIRST_UPI, -- UPI
             CURRENT_TIMESTAMP, -- timestamp
             USER,  -- user stamp
             v_crc(1), -- crc64
             SEQ_LENGTH, -- length
             v_seq(1), -- sequence
             NULL, -- long sequence
             v_md5(1) -- md5
             );
    COMMIT;
  END initialize_rna_table;

  /* Insert the first record. */
  PROCEDURE initialize_xref_table AS
  BEGIN
    INSERT
      INTO RNACEN.xref
      VALUES(1,                 -- id
             1,                 -- dbid
             1,                 -- first release
             FIRST_UPI,         -- UPI
             1,                 -- last release
             'N',               -- deleted
             CURRENT_TIMESTAMP, -- timestamp
             USER,              -- userstamp
             v_acc(1),          -- accession
             1,                 -- version_I
             'gi',              -- GI
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
    -- add first release
    INSERT
      INTO RNACEN.rnc_release
      VALUES(1,                 -- id
             1,                 -- dbid
             CURRENT_TIMESTAMP, -- release date
             'F',               -- release type
             'L',               -- release status
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
    Import new data into the staging table
  */
  PROCEDURE import_staging_data(
    p_test_id NUMBER
  ) AS
  BEGIN
    INSERT INTO
      RNACEN.load_rnacentral_test
      VALUES(v_crc(p_test_id), -- crc64
             SEQ_LENGTH,       -- length
             v_seq(p_test_id), -- seq
             NULL,             -- long seq
             v_acc(p_test_id), -- accession
             v_ver(p_test_id), -- version
             v_tax(p_test_id), -- taxid
             v_md5(p_test_id)  --md5
             );
    COMMIT;
  END import_staging_data;

  /*
    Truncate all test tables
  */
  PROCEDURE teardown AS
  BEGIN

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

    EXECUTE IMMEDIATE 'ALTER TABLE RNACEN.rnc_release DISABLE CONSTRAINT ID_PK';

    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rnc_release DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.rna DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.xref DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE RNACEN.load_rnacentral_test DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_DELETED DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_NOT_DELETED DROP STORAGE';

    EXECUTE IMMEDIATE 'ALTER TABLE RNACEN.rnc_release ENABLE CONSTRAINT ID_PK';

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
    Inspired by https://github.com/uzilan/asserts-package
  */
  PROCEDURE assertEquals(
    proc IN VARCHAR2,
    expected in NUMBER,
    actual in NUMBER)
  IS
  BEGIN
    IF NOT expected = actual THEN
      DBMS_OUTPUT.put_line(proc || 'expected ' || expected || ', got ' || actual);
    ELSE
      DBMS_OUTPUT.put_line(proc || ' ok');
    END IF;
  end assertEquals;

  /*
    Same sequence.
  */
  PROCEDURE check_test1 AS
    l_count NUMBER;
  BEGIN

    SELECT count(*) INTO l_count FROM rna;
    assertEquals('check_test1', 1, l_count);

    SELECT count(*) INTO l_count FROM xref;
    assertEquals('check_test1', 2, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'Y' AND created = 1 AND LAST = 1;
    assertEquals('check_test1', 1, l_count);

    SELECT count(*) INTO l_count FROM xref
    WHERE upi = FIRST_UPI AND deleted = 'N' AND created = 1 AND LAST = 2;
    assertEquals('check_test1', 1, l_count);

  END check_test1;

  /*
    Check results
  */
  PROCEDURE check_result(
    p_test_id NUMBER
  ) AS
  BEGIN

    -- todo: formalize this better
    CASE p_test_id
      WHEN 2 THEN check_test1;
      ELSE
        NULL;
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
    import_staging_data( p_test_id );

    -- run update
    LOAD_RETRO_SETBASED.load_job_setbased;

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

    -- truncate tables in the beginning
    teardown;

    FOR l_cntr IN 2..17
    LOOP
      DBMS_OUTPUT.put_line('Running test ' || l_cntr);
      run_test(l_cntr);
    END LOOP;

  END run_tests;


END RNC_TEST;