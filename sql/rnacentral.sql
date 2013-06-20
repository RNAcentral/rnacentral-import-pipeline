CREATE OR REPLACE
PACKAGE BODY "LOAD_RETRO_SETBASED"
IS

/*
* Loads the data from the staging table LOAD_RNACENTRAL_TEST
* into a temporary table LOAD_RETRO_TMP.
* Retrieves the UPI of the corresponding RNA
* when the new and existing sequences are the same.
*/
PROCEDURE load_retro_tmp_table(
    p_in_dbid         IN RNACEN.RNC_DATABASE.ID%TYPE,
    p_in_load_release IN RNACEN.RNC_RELEASE.ID%TYPE )
IS
BEGIN

  BEGIN
    EXECUTE IMMEDIATE 'DROP INDEX LOAD_RETRO_TMP$AC_DBID_UPI';
  EXCEPTION
    WHEN OTHERS THEN
      NULL;
  END;

  EXECUTE IMMEDIATE 'TRUNCATE TABLE LOAD_RETRO_TMP DROP STORAGE';
  EXECUTE IMMEDIATE 'ALTER SESSION FORCE PARALLEL DML';

  EXECUTE IMMEDIATE '
INSERT
  /*+ APPEND PARALLEL */
INTO
  LOAD_RETRO_TMP
  (
    IN_DBID ,
    IN_LOAD_RELEASE ,
    IN_CRC64 ,
    IN_LEN ,
    IN_SEQ_SHORT ,
    IN_SEQ_LONG ,
    IN_AC ,
    IN_VERSION ,
    IN_MD5 ,
    IN_TAXID ,
    COMPARABLE_PROT_UPI
  )
SELECT
  /*+ PARALLEL */
  :p_in_dbid,
  :p_in_load_release,
  in_crc64,
  in_len,
  in_seq_short,
  in_seq_long,
  in_ac,
  in_version,
  IN_MD5,
  in_taxid,
  CASE
    WHEN prot_upi IS NOT NULL
    THEN -- a rna with the same md5 exists
      CASE
        WHEN
          (
            (
              in_len   > 4000
            AND in_len = prot_len
            )
          AND dbms_lob.compare( in_seq_long, prot_seq_long) = 0
          )
        THEN -- a rna with the same md5, len, sequence (lob) exists
          prot_upi
        WHEN
          (
            (
              in_len  <= 4000
            AND in_len = prot_len
            )
          AND in_seq_short = prot_seq_short
          )
        THEN -- a rna with the same md5, len, sequence (not lob) exists
          prot_upi
        ELSE -- the rna is not really the same
          NULL
      END
  END comparable_prot_upi
FROM
  (
  WITH
    -- `distinct_loaded_rows` contains sequences from the staging table
    distinct_loaded_rows AS
    (
      SELECT DISTINCT
        CRC64,
        LEN,
        SEQ_SHORT,
        useful_clob_object(seq_long) MY_SEQ_LONG,
        AC,
        VERSION,
        TAXID,
        MD5
      FROM
        RNACEN.LOAD_RNACENTRAL_TEST
     )
  SELECT
    l.crc64 in_crc64,
    L.LEN in_len,
    L.SEQ_SHORT in_seq_short,
    TREAT(l.MY_SEQ_LONG AS USEFUL_CLOB_OBJECT).UCO AS in_seq_long,
    l.ac in_ac,
    l.version in_version,
    L.MD5 IN_MD5,
    l.taxid in_taxid,
    p.id prot_id,
    p.seq_short prot_seq_short,
    p.seq_long prot_seq_long,
    p.len prot_len,
    P.UPI PROT_UPI
  FROM
    distinct_loaded_rows l,
    RNACEN.rna p
  WHERE p.md5 (+) = l.md5)' USING TO_CHAR (
    p_in_dbid), TO_CHAR (p_in_load_release);

  COMMIT;

  BEGIN
    DBMS_STATS.GATHER_TABLE_STATS ( ownname => 'RNACEN', tabname =>
    'LOAD_RETRO_TMP', estimate_percent => 1);
  END;

END load_retro_tmp_table;

/*
*
*
*
*/
PROCEDURE load_md5_stats_table IS
BEGIN

  EXECUTE IMMEDIATE 'TRUNCATE TABLE load_md5_stats DROP STORAGE';
  INSERT
    /*+ APPEND PARALLEL (LOAD_MD5_STATS 4) */
  INTO
    RNACEN.LOAD_MD5_STATS
    (
      IN_MD5,
      CNT,
      CNT_DST_SEQ_SHORT,
      CNT_DST_SEQ_LONG
    )
  SELECT
    IN_MD5,
    COUNT (*) cnt,
    COUNT (DISTINCT IN_SEQ_SHORT) CNT_DST_SEQ_SHORT,
    COUNT (DISTINCT DBMS_LOB.SUBSTR (IN_SEQ_LONG, 4000, 1)) CNT_DST_SEQ_LONG
  FROM
    (
      SELECT
        /*+ PARALLEL (L 4) */
        l.in_seq_short,
        l.in_seq_long,
        l.in_md5
      FROM
        LOAD_RETRO_TMP L
      WHERE
        l.comparable_prot_upi IS NULL
    )
  GROUP BY
    IN_MD5;

  COMMIT;

END load_md5_stats_table;

/*
* This table is never actually used in the algorithm.
* Maintained for logging purposes only.
*/
PROCEDURE load_md5_collisions_table IS
BEGIN

  EXECUTE IMMEDIATE 'TRUNCATE TABLE load_md5_collisions DROP STORAGE';

  INSERT
    /*+ APPEND PARALLEL (LOAD_MD5_COLLISIONS 4) */
  INTO
    RNACEN.LOAD_MD5_COLLISIONS
    (
      IN_MD5,
      CNT,
      CNT_DST_SEQ_SHORT,
      CNT_DST_SEQ_LONG
    )
  SELECT
    /*+ PARALLEL (LOAD_MD5_STATS 4) */
    IN_MD5,
    cnt,
    CNT_DST_SEQ_SHORT,
    CNT_DST_SEQ_LONG
  FROM
    LOAD_MD5_STATS
  WHERE
    (
      CNT_DST_SEQ_SHORT > 1
    OR CNT_DST_SEQ_LONG > 1
    );

  COMMIT;

END load_md5_collisions_table;

/*
* Loads previously unobserved sequences into the temporary table.
* Assigns new ids by calling Java via a PL/SQL wrapper package.
* SEQ_UPI is the sequence object.
*/
PROCEDURE load_md5_new_sequences_table IS
BEGIN

  EXECUTE IMMEDIATE 'TRUNCATE TABLE load_md5_new_sequences DROP STORAGE';

  INSERT
    /*+ APPEND PARALLEL (load_md5_new_sequences 4) */
  INTO
    RNACEN.load_md5_new_sequences
    (
      IN_MD5,
      PROT_ID,
      PROT_UPI
    )
  SELECT
    /*+ PARALLEL (LOAD_MD5_STATS 4) */
    IN_MD5,
    SEQ_UPI.NEXTVAL PROT_ID,
    CAST (RNACEN.UPI.GET_UPI(SEQ_UPI.CURRVAL) AS CHAR(13)) PROT_UPI
  FROM
    LOAD_MD5_STATS
  WHERE
    NOT
    (
      CNT_DST_SEQ_SHORT > 1 OR CNT_DST_SEQ_LONG > 1
    );
  COMMIT;

  DBMS_STATS.GATHER_TABLE_STATS ( ownname => 'RNACEN', tabname =>
  'load_md5_new_sequences', estimate_percent => 10 );

END load_md5_new_sequences_table;

/*
* Propagate the newly assigned identifiers
* to the field COMPARABLE_PROT_UPI.
*/
PROCEDURE set_comparable_prot_upi IS
BEGIN

  -- Now COMPARABLE_PROT_UPI IS NOT NULL also for new sequences

  UPDATE
    /*+ PARALLEL (L 4) */
    LOAD_RETRO_TMP L
  SET
    L.COMPARABLE_PROT_UPI =
    (
      SELECT
        /*+ PARALLEL (N 4) */
        N.PROT_UPI
      FROM
        load_md5_new_sequences N
      WHERE
        N.IN_MD5 = L.IN_MD5
    )
  WHERE
    L.COMPARABLE_PROT_UPI IS NULL;

  COMMIT;

  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."LOAD_RETRO_TMP$AC_DBID_UPI" ON LOAD_RETRO_TMP (IN_AC, IN_DBID, COMPARABLE_PROT_UPI) PARALLEL TABLESPACE "RNACEN_IND"';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."LOAD_RETRO_TMP$AC_DBID_UPI" NOPARALLEL';
  BEGIN
    DBMS_STATS.GATHER_TABLE_STATS ( ownname => 'RNACEN', tabname =>
    'LOAD_RETRO_TMP', estimate_percent => 10 );
  END;

END;

/*
* Deposit the new sequences into the main table
* containing all sequences from all sources.
*/
PROCEDURE store_new_sequences IS
BEGIN

  INSERT
    /*+ APPEND PARALLEL (RNA 2) */
  INTO
    RNACEN.rna
    (
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
    PROT_ID,
    PROT_UPI,
    in_crc64,
    in_len,
    in_seq_short,
    IN_SEQ_LONG,
    IN_MD5,
    SYSDATE,
    USER
  FROM
    -- Using ROW_NUMBER because seq_long is a LOB
    -- can not use DISTINCT with LOBs

    -- ROW_NUMBER is an analytic function
    (
      SELECT
        /*+ PARALLEL (L 2) PARALLEL (N 2) */
        n.PROT_ID,
        n.PROT_UPI,
        l.in_crc64,
        l.in_len,
        l.in_seq_short,
        l.IN_SEQ_LONG,
        l.in_md5,
        ROW_NUMBER () OVER (PARTITION BY n.prot_upi ORDER BY NULL) rn
      FROM
        load_retro_tmp l,
        load_md5_new_sequences n
      WHERE
        n.in_md5 = l.in_md5
        -- AND l.comparable_prot_upi IS NULL
    )
  WHERE
    rn <= 1;

  COMMIT;

END;

-------------------
-- XREF analysis --
-------------------

/*
*
*
*
*/
PROCEDURE prepare_xref_pel_tables IS
BEGIN

    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_DELETED DROP STORAGE';
    EXECUTE IMMEDIATE 'TRUNCATE TABLE XREF_PEL_NOT_DELETED DROP STORAGE';
    BEGIN
      -- INCAPSULATED TO IGNORE ANY 'DOES NOT EXISTS ERROR'
      -- WHEN RESTARTING THE PROCEDURE AFTER UNPREDICTABLE EXCEPTIONS
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$AC
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$AC"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$CREATED
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$CREATED"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$GI
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$GI"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$UPI
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$UPI"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$TAXID
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$TAXID"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_DELETED$AC$UPPER
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_DELETED$AC$UPPER"';
      --------------------------------------------------------
      --  Constraints for Table XREF
      --------------------------------------------------------
      /*
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "
      XREF_PEL_DELETED_CK1";
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DBID" NOT NULL DISABLE
      NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("CREATED" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("LAST" NOT NULL DISABLE
      NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("UPI" NOT NULL DISABLE
      NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("VERSION_I" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DELETED" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("TIMESTAMP" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("USERSTAMP" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("AC" NOT NULL DISABLE
      NOVALIDATE);
      */
      --------------------------------------------------------
      --  Ref Constraints for Table XREF
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK1"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK2"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK3"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" DROP CONSTRAINT "XREF_PEL_DELETED_FK4"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$AC
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$CREATED
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$CREATED"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$GI
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$GI"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$UPI
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$UPI"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$TAXID
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$TAXID"';
      --------------------------------------------------------
      --  DDL for Index XREF_PEL_NOT_DELETED$AC$UPPER
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'DROP INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC$UPPER"';
      --------------------------------------------------------
      --  Constraints for Table XREF
      --------------------------------------------------------
      /*
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "
      XREF_PEL_NOT_DELETED_CK1";
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DBID" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("CREATED" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("LAST" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("UPI" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("VERSION_I" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DELETED" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("TIMESTAMP" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("USERSTAMP" NOT NULL
      DISABLE NOVALIDATE);
      ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("AC" NOT NULL
      DISABLE NOVALIDATE);
      */
      --------------------------------------------------------
      --  Ref Constraints for Table XREF
      --------------------------------------------------------
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK1"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK2"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK3"';
      EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" DROP CONSTRAINT "XREF_PEL_NOT_DELETED_FK4"';

    EXCEPTION
    WHEN OTHERS THEN
      NULL;
    END;

END prepare_xref_pel_tables;

/*
    -- affects existing xrefs with unchanged UPI / VERSION
    -- OR cancels old entries (UPI / CHANGED VERSION)
    -- new entries are inserted in the query below

    -- XREF (AC,DBID) ESISTENTI
    -- DA TOCCARE SOLTANTO (UPI/VERSION INALTERATI)
    -- O DI CUI SI CANCELLA LA VECCHIA IMMAGINE (UPI/VERSION CAMBIATI)
    -- (INSERIMENTO NUOVA IMMAGINE IN QUERY SOTTOSTANTE)
*
* Copy over existing active xrefs to the PEL tables.
* Retire entries which have different UPIs or updated versions.
* Looks only at existing entries that match something in the staging table
*/
PROCEDURE populate_pel_tables1 (
  v_previous_release RNACEN.rnc_release.id%TYPE)
IS
BEGIN

  INSERT --+ APPEND PARALLEL (XREF_PEL_DELETED 2) PARALLEL (
    -- XREF_PEL_NOT_DELETED 2)
    ALL
    WHEN DELETED = 'Y' THEN
  INTO
    XREF_PEL_DELETED
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
    VALUES
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
    WHEN DELETED = 'N' THEN
  INTO
    XREF_PEL_NOT_DELETED
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
    VALUES
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
  SELECT
    /*+ PARALLEL (X 2) PARALLEL (L 2) */
    x.DBID,
    x.CREATED,
    x.UPI,
    x.VERSION_I,
    x."TIMESTAMP",
    x.USERSTAMP,
    x.AC,
    x."VERSION",
    x.GI,
    -- x.LAST
    CASE
      WHEN
        (
          x.LAST  < l.in_load_release -- Question: isn't this always true?
        AND X.UPI = L.COMPARABLE_PROT_UPI
        AND
          (
            X."VERSION" = L.IN_VERSION -- the external db version matches what's already in xref
          OR
            (
              X."VERSION"    IS NULL
            AND L.IN_VERSION IS NULL
            )
          )
        )
      THEN L.IN_LOAD_RELEASE
      ELSE NVL (v_previous_release, x.LAST)
    END "LAST",
    -- x.DELETED
    CASE
      WHEN
        (
          x.LAST  < l.in_load_release
        AND X.UPI = L.COMPARABLE_PROT_UPI
        AND
          (
            X."VERSION" = L.IN_VERSION
          OR
            (
              X."VERSION"    IS NULL
            AND L.IN_VERSION IS NULL
            )
          )
        )
      THEN 'N'
      ELSE 'Y'
    END "DELETED",
    -- x.taxid
    CASE
      WHEN
        (
          x.LAST  < l.in_load_release
        AND X.UPI = L.COMPARABLE_PROT_UPI
        AND
          (
            X."VERSION" = L.IN_VERSION
          OR
            (
              X."VERSION"    IS NULL
            AND L.IN_VERSION IS NULL
            )
          )
        )
      THEN NVL (L.IN_TAXID, x.TAXID) -- override taxid
      ELSE x.TAXID
    END "TAXID"
  FROM
    LOAD_RETRO_TMP l,
    XREF x
  WHERE
      x.AC                   = l.IN_AC
  AND x.DBID                 = l.IN_DBID
  AND l.comparable_prot_upi IS NOT NULL
  AND x.deleted              = 'N'; -- select only active existing xrefs

  COMMIT;

END populate_pel_tables1;

/*
* it looks like this query inserts updated xrefs
* which replace xrefs retired in populate_pel_tables1

* deals only with XREF_PEL_NOT_DELETED
*/

PROCEDURE populate_pel_tables2 IS
BEGIN

-- affects existing XREFs
-- for which new xrefs are inserted (UPI / VERSION CHANGED)
-- (previous xrefs were canceled in the previous query)

-- XREF (AC,DBID) ESISTENTI
-- DI CUI SI INSERISCE LA NUOVA IMMAGINE (UPI/VERSION CAMBIATI)
-- (CANCELLAZIONE IMMAGINE PRECEDENTE IN QUERY SOPRASTANTE)

-- adds new versions of the existing entries


  INSERT --+ APPEND PARALLEL (XREF_PEL_NOT_DELETED 8)
  INTO
    XREF_PEL_NOT_DELETED
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
  SELECT
    /*+ PARALLEL (X 4) PARALLEL (L 4) */
    l.IN_DBID,
    l.IN_LOAD_RELEASE CREATED,
    l.COMPARABLE_PROT_UPI UPI,
    CASE
      WHEN
        (
          x.UPI != l.COMPARABLE_PROT_UPI
        )
      THEN x.VERSION_I + 1
      ELSE X.VERSION_I
    END VERSION_I,
    SYSDATE "TIMESTAMP",
    USER USERSTAMP,
    L.IN_AC AC,
    l.in_version "VERSION",
    x.GI,
    l.IN_LOAD_RELEASE "LAST",
    'N' DELETED,                       -- Question: unusual syntax
    IN_TAXID TAXID
  FROM
    LOAD_RETRO_TMP l,
    XREF x
  WHERE
      x.AC                   = L.IN_AC
  AND x.DBID                 = L.IN_DBID
  AND l.COMPARABLE_PROT_UPI IS NOT NULL
  -- the next condition differentiates this procedure from populate_pel_tables1
  AND NOT
    (
      x.LAST  < l.IN_LOAD_RELEASE
    AND x.UPI = l.COMPARABLE_PROT_UPI
    AND
      (
        X."VERSION" = L.IN_VERSION
      OR
        (
          X."VERSION"    IS NULL
        AND L.IN_VERSION IS NULL
        )
      )
    )
  AND x.deleted = 'N';

  COMMIT;

END populate_pel_tables2;

/*
* Retrieve max versions for all xrefs in the staging table, which are new
* (assigned version = 0) or which are inactive because have been deleted,
* but should be activated now. UPIs are kept.
*/
PROCEDURE load_upi_max_versions_table(
  p_in_dbid IN RNACEN.RNC_DATABASE.ID%TYPE
) IS
BEGIN

  EXECUTE IMMEDIATE 'TRUNCATE TABLE load_upi_max_versions DROP STORAGE';

  -- inserts a new xref
  -- for dbid, ac where xref:
  -- a) does not exists at all (VERSION_I = 0 will become 1 when inserting
  --    notice outer join to previous_xref
  -- b) does not exists because has been deleted
  --    1) and was having the same UPI, Version (VERSION_I = OLD_VERSION_I)
  --    2) and was having a different UPI, VERSION (VERSION_I = OLD_VERSION_I + 1)
  INSERT /*+ APPEND */ INTO load_upi_max_versions
  SELECT
    /*+ PARALLEL (PREVIOUS_XREF 4) PARALLEL (L 4) */
    DISTINCT L.IN_AC,
    L.IN_DBID,
    MAX(NVL (PREVIOUS_XREF.VERSION_I, 0)) MAX_VERSION_I,
    PREVIOUS_XREF.UPI
  FROM
    LOAD_RETRO_TMP L,
    RNACEN.XREF PREVIOUS_XREF
  WHERE
    L.COMPARABLE_PROT_UPI IS NOT NULL  -- Question: is it not always null?
  AND L.IN_DBID = p_in_dbid
  AND PREVIOUS_XREF.DBID (+)  = L.IN_DBID
  AND PREVIOUS_XREF.AC   (+)  = L.IN_AC
  AND NOT EXISTS -- no active entry exists
    (
      SELECT
        /*+ PARALLEL (X 4) */
        1
      FROM
        RNACEN.XREF X
      WHERE
          X.AC      = L.IN_AC
      AND X.DBID    = L.IN_DBID
      AND X.DELETED = 'N'
    )
  GROUP BY
    L.IN_AC,
    L.IN_DBID,
    PREVIOUS_XREF.UPI;

  COMMIT;

  BEGIN
    DBMS_STATS.GATHER_TABLE_STATS ( OWNNAME => 'RNACEN', TABNAME => 'load_upi_max_versions', ESTIMATE_PERCENT => 10 );
  END;

END load_upi_max_versions_table;

/*
* Get absolute maximum version for each accession number/db id pair
*/
PROCEDURE load_max_versions_table IS
BEGIN
  EXECUTE IMMEDIATE 'TRUNCATE TABLE load_max_versions DROP STORAGE';

  INSERT /*+ APPEND */ INTO load_max_versions
  SELECT
    /*+ PARALLEL (L 4) */
    DISTINCT l.AC,
    L.DBID,
    MAX(L.MAX_VERSION_I) MAX_VERSION_I
  FROM load_upi_max_versions l
  GROUP BY
    L.AC,
    L.DBID;

  COMMIT;

  BEGIN
    DBMS_STATS.GATHER_TABLE_STATS ( OWNNAME => 'RNACEN', TABNAME => 'load_max_versions', ESTIMATE_PERCENT => 10 );
  END;

END load_max_versions_table;

/*
*
*
*
*/
PROCEDURE populate_pel_tables3 (
  p_in_dbid IN RNACEN.RNC_DATABASE.ID%TYPE
) IS
BEGIN

  INSERT --+ APPEND PARALLEL (XREF_PEL_NOT_DELETED 4)
  INTO
    RNACEN.XREF_PEL_NOT_DELETED
    (
      AC,
      DBID,
      "VERSION",
      VERSION_I,
      UPI,
      CREATED,
      "LAST",
      "DELETED",
      TAXID,
      GI,
      "TIMESTAMP",
      USERSTAMP
    )
  SELECT
    /*+ PARALLEL (T 4) PARALLEL (LUMV 4) */
    T.IN_AC,
    T.IN_DBID,
    T.IN_VERSION,
    CASE
      WHEN T.MAX_PREVIOUS_XREF_VERSION_I = 0
      THEN 1                                  -- 1 for new xrefs
      WHEN LUMV.UPI = T.COMPARABLE_PROT_UPI
      THEN T.MAX_PREVIOUS_XREF_VERSION_I      -- previous max, same xref is reactivated
      ELSE T.MAX_PREVIOUS_XREF_VERSION_I + 1  -- bump up the version
    END VERSION_I,
    T.COMPARABLE_PROT_UPI UPI,
    T.IN_LOAD_RELEASE CREATED,
    T.IN_LOAD_RELEASE LAST,
    'N' DELETED,
    T.IN_TAXID,
    NULL GI,
    SYSDATE "TIMESTAMP",
    USER USERSTAMP
  FROM
    LOAD_UPI_MAX_VERSIONS lumv,
    (
      SELECT
        /*+ PARALLEL (L 4) PARALLEL (LMV 4) */
        L.IN_AC,
        L.IN_DBID,
        L.IN_VERSION,
        LMV.MAX_VERSION_I MAX_PREVIOUS_XREF_VERSION_I,
        L.COMPARABLE_PROT_UPI,
        L.IN_LOAD_RELEASE,
        L.IN_TAXID
      FROM
        LOAD_RETRO_TMP L,
        LOAD_MAX_VERSIONS LMV
      WHERE
        L.COMPARABLE_PROT_UPI IS NOT NULL
      AND L.IN_DBID = p_in_dbid
      AND LMV.AC           = L.IN_AC
      AND LMV.DBID         = L.IN_DBID
    )
    T
  WHERE
      LUMV.AC (+)            = T.IN_AC
  AND LUMV.MAX_VERSION_I (+) = T.MAX_PREVIOUS_XREF_VERSION_I
  AND LUMV.DBID (+)          = T.IN_DBID;

  COMMIT;

END populate_pel_tables3;

/*
* deals only with XREF_PEL_DELETED
*
*
*/
PROCEDURE populate_pel_tables4 (
  p_in_dbid IN RNACEN.RNC_DATABASE.ID%TYPE,
  v_previous_release IN RNACEN.RNC_RELEASE.ID%TYPE
) IS
BEGIN

  INSERT --+ APPEND PARALLEL (XREF_PEL_DELETED 4)
  INTO
    XREF_PEL_DELETED
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
  SELECT
    /*+ PARALLEL (X 4) */
    -- XREF (AC,DBID) ESISTENTI
    -- NON PRESENTI NELLA RELEASE
    --       PER I QUALI DOBBIAMO CANCELLARE LA VECCHIA IMMAGINE
    -- O CON IMMAGINI GIA' CANCELLATE
    --       DA MANTENERE INALTERATE
    x.DBID,
    x.CREATED,
    x.UPI,
    x.VERSION_I,
    x.TIMESTAMP,
    x.USERSTAMP,
    x.AC,
    x.VERSION,
    x.gi,
    CASE X.DELETED
      WHEN 'N'
      THEN
        -- IMMAGINE DA CANCELLARE (RELEASE PRECEDENTE ULTIMA VALIDA)
        NVL (v_previous_release, X.LAST)
      ELSE -- IMMAGINE GIA' CANCELLATA RIMANE INALTERATA
        X.LAST
    END LAST,
    'Y' DELETED,
    x.taxid
  FROM
    RNACEN.XREF X
  WHERE
    X.dbid = p_in_dbid
  AND
    (
      x.deleted = 'Y'
    OR NOT EXISTS
      (
        SELECT
          /*+ PARALLEL (L 4) */
          1
        FROM
          LOAD_RETRO_TMP L
        WHERE
          L.COMPARABLE_PROT_UPI IS NOT NULL
        AND X.AC                 = L.IN_AC
        AND X.DBID               = L.IN_DBID
      )
    );

  COMMIT;

END populate_pel_tables4;

/*
* Copy over entries that had been deleted in previous releases
* or that don't match anything in the staging table.
* The entries must be copied because the old xref table is
* exchanged with the new XREF_PEL tables.
*/
PROCEDURE populate_pel_tables5 (
  p_in_dbid IN RNACEN.RNC_DATABASE.ID%TYPE,
  v_previous_release IN RNACEN.RNC_RELEASE.ID%TYPE,
  v_release_type IN RNACEN.RNC_RELEASE.release_type%TYPE
) IS
BEGIN

  INSERT --+ APPEND PARALLEL (XREF_PEL_DELETED 2) PARALLEL (
    -- XREF_PEL_NOT_DELETED 2)
    ALL
    WHEN DELETED = 'Y' THEN
  INTO
    XREF_PEL_DELETED
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
    VALUES
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
    WHEN DELETED = 'N' THEN
  INTO
    XREF_PEL_NOT_DELETED
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
    VALUES
    (
      DBID,
      CREATED,
      UPI,
      VERSION_I,
      "TIMESTAMP",
      USERSTAMP,
      AC,
      "VERSION",
      GI,
      "LAST",
      "DELETED",
      "TAXID"
    )
  SELECT
    /*+ PARALLEL (X 4) */
    -- XREF (AC,DBID) ESISTENTI
    -- NON PRESENTI NELLA RELEASE PER I QUALI
    --       DOBBIAMO MANTENERE LA VECCHIA IMMAGINE INALTERATA
    --       (nel caso di release INCREMENTALE)
    --     OPPURE
    --       CANCELLARE LA VECCHIA IMMAGINE
    --       (nel caso di release CUMULATIVA)
    -- O CON IMMAGINI GIA' CANCELLATE
    --       DA MANTENERE INALTERATE
    x.DBID,
    x.CREATED,
    x.UPI,
    x.VERSION_I,
    x.TIMESTAMP,
    x.USERSTAMP,
    x.AC,
    x.VERSION,
    x.gi,
    x.LAST,
    CASE X.DELETED
      WHEN 'N'
      THEN
        CASE
            -- A cumulative release is substituted completely
            -- by a IMMEDIATE successive cumulative release
          WHEN
            (
              x.LAST                    = v_previous_release
            AND v_previous_release     IS NOT NULL
            AND v_release_type          = 'C'
            --AND v_previous_release_type = 'C'
            )
          THEN
            -- IMMAGINE DA CANCELLARE (RELEASE PRECEDENTE ULTIMA VALIDA)
            'Y'
          ELSE
            -- IMMAGINE RIMANE INALTERATA
            X.DELETED
        END
      ELSE
        -- IMMAGINE GIA' CANCELLATA RIMANE INALTERATA
        X.DELETED
    END DELETED,
    x.taxid
  FROM
    RNACEN.XREF X
  WHERE
    X.DBID = p_in_dbid
  AND
    (
      X.DELETED = 'Y'
    OR NOT EXISTS -- the entries in xref don't match anything in the staging table
      (
        SELECT
          /*+ PARALLEL (L 4) */
          1
        FROM
          LOAD_RETRO_TMP L
        WHERE
          L.COMPARABLE_PROT_UPI IS NOT NULL
        AND X.AC                 = L.IN_AC
        AND X.DBID               = L.IN_DBID
      )
    );

  COMMIT;

END populate_pel_tables5;


/*
*
*
*
*/
PROCEDURE prepare_pel_exchange IS
BEGIN

  --------------------------------------------------------
  --  DDL for Index XREF_PEL_DELETED$AC
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$AC" ON "RNACEN"."XREF_PEL_DELETED" ("AC") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$AC" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_DELETED$CREATED
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$CREATED" ON "RNACEN"."XREF_PEL_DELETED" ("CREATED") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$CREATED" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_DELETED$GI
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$GI" ON "RNACEN"."XREF_PEL_DELETED" ("GI") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$GI" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_DELETED$UPI
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$UPI" ON "RNACEN"."XREF_PEL_DELETED" ("UPI") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$UPI" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_DELETED$TAXID
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$TAXID" ON "RNACEN"."XREF_PEL_DELETED" ("TAXID") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$TAXID" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_DELETED$AC$UPPER
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_DELETED$AC$UPPER" ON "RNACEN"."XREF_PEL_DELETED" (UPPER("AC")) PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_DELETED$AC$UPPER" NOPARALLEL';
  --------------------------------------------------------
  --  Constraints for Table XREF
  --------------------------------------------------------
  /*
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_CK1" CHECK (deleted IN ('Y', 'N')) ENABLE VALIDATE;
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DBID" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("CREATED" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("LAST" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("UPI" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("VERSION_I" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("DELETED" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("TIMESTAMP" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("USERSTAMP" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_DELETED" MODIFY ("AC" NOT NULL ENABLE VALIDATE);
  */
  --------------------------------------------------------
  --  Ref Constraints for Table XREF
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK1" FOREIGN KEY ("CREATED") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK2" FOREIGN KEY ("DBID") REFERENCES "RNACEN"."RNC_DATABASE" ("ID") ENABLE VALIDATE';

  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK3" FOREIGN KEY ("LAST") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_DELETED" ADD CONSTRAINT "XREF_PEL_DELETED_FK4" FOREIGN KEY ("UPI") REFERENCES "RNACEN"."RNA" ("UPI") ENABLE VALIDATE';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_NOT_DELETED$AC
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("AC") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_NOT_DELETED$CREATED
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$CREATED" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("CREATED") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$CREATED" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_NOT_DELETED$GI
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$GI" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("GI") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$GI" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_NOT_DELETED$UPI
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$UPI" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("UPI") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$UPI" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_NOT_DELETED$TAXID
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$TAXID" ON "RNACEN"."XREF_PEL_NOT_DELETED" ("TAXID") PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$TAXID" NOPARALLEL';
  --------------------------------------------------------
  --  DDL for Index XREF_PEL_NOT_DELETED$AC$UPPER
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'CREATE INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC$UPPER" ON "RNACEN"."XREF_PEL_NOT_DELETED" (UPPER("AC")) PARALLEL TABLESPACE "RNACEN_IND" COMPRESS 1';
  EXECUTE IMMEDIATE 'ALTER INDEX "RNACEN"."XREF_PEL_NOT_DELETED$AC$UPPER" NOPARALLEL';
  --------------------------------------------------------
  --  Constraints for Table XREF
  --------------------------------------------------------
  /*
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_CK1" CHECK (deleted IN ('Y', 'N')) ENABLE VALIDATE;
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DBID" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("CREATED" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("LAST" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("UPI" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("VERSION_I" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("DELETED" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("TIMESTAMP" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("USERSTAMP" NOT NULL ENABLE VALIDATE);
  ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" MODIFY ("AC" NOT NULL ENABLE VALIDATE);
  */
  --------------------------------------------------------
  --  Ref Constraints for Table XREF
  --------------------------------------------------------
  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK1" FOREIGN KEY ("CREATED") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK2" FOREIGN KEY ("DBID") REFERENCES "RNACEN"."RNC_DATABASE" ("ID") ENABLE VALIDATE';
  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK3" FOREIGN KEY ("LAST") REFERENCES "RNACEN"."RNC_RELEASE" ("ID") ENABLE VALIDATE';
  EXECUTE IMMEDIATE 'ALTER TABLE "RNACEN"."XREF_PEL_NOT_DELETED" ADD CONSTRAINT "XREF_PEL_NOT_DELETED_FK4" FOREIGN KEY ("UPI") REFERENCES "RNACEN"."RNA" ("UPI") ENABLE VALIDATE';

  -- THIS ELIMINATES THE NEED OF CALCULATING THE CBO STATISTICS ON THE
  -- SUBPARTITION AFTER THE EXCHANGE
  DBMS_STATS.GATHER_TABLE_STATS (OWNNAME => 'RNACEN', TABNAME =>
  'XREF_PEL_DELETED', ESTIMATE_PERCENT => DBMS_STATS.AUTO_SAMPLE_SIZE,
  METHOD_OPT => 'FOR ALL COLUMNS SIZE AUTO', DEGREE => 8, CASCADE => TRUE);

  -- THIS ELIMINATES THE NEED OF CALCULATING THE CBO STATISTICS ON THE
  -- SUBPARTITION AFTER THE EXCHANGE
  DBMS_STATS.GATHER_TABLE_STATS (OWNNAME => 'RNACEN', TABNAME =>
  'XREF_PEL_NOT_DELETED', ESTIMATE_PERCENT => DBMS_STATS.AUTO_SAMPLE_SIZE,
  METHOD_OPT => 'FOR ALL COLUMNS SIZE AUTO', DEGREE => 4, CASCADE => TRUE);

END prepare_pel_exchange;




PROCEDURE incremental_release (
  v_previous_release IN RNACEN.rnc_release.id%TYPE,
  p_in_dbid          IN RNACEN.RNC_DATABASE.ID%TYPE
)
IS
BEGIN

  /*
  -- this is an INCREMENTAL release
  -- this is NOT a FULL release
  -- better to AVOID Partition Exchange Loading
  -- updates last, deleted and taxid
  -- of existing xrefs

  Unfortunately when using parallel with conventional DML we are hitting a
  bug. I am disabling parallel DML but leaving PARALLEL hints in place so that
  the query part of the DML statements should still be able to use it.
  Bug 9831227 : Parallel DML fails with ORA-12839 or ORA-600
  Description
  If you have a parent table and child table with a parent referential
  constraint
  then running Parallel DML again this may fail with an unexpected ORA-12839
  or even fail with an ORA-600.
  Workaround
  Setting "_disable_parallel_conventional_load" = true can help avoid this.
  */

  EXECUTE IMMEDIATE 'ALTER SESSION DISABLE PARALLEL DML';
  /*
  -- updates last, deleted and taxid
  -- of existing xrefs
  */
  UPDATE
    /*+ PARALLEL (U 4) */
    RNACEN.XREF U
  SET
    (
      U.LAST,
      U.DELETED,
      u.taxid
    )
    =
    (
      SELECT
        /*+ PARALLEL (L 4) */
        CASE
          WHEN
            (
              U.UPI = L.COMPARABLE_PROT_UPI
            AND
              (
                U."VERSION" = L.IN_VERSION
              OR
                (
                  U."VERSION"    IS NULL
                AND L.IN_VERSION IS NULL
                )
              )
            )
          THEN L.IN_LOAD_RELEASE
          ELSE v_previous_release
        END "LAST",
        CASE
          WHEN
            (
              U.UPI = L.COMPARABLE_PROT_UPI
            AND
              (
                U."VERSION" = L.IN_VERSION
              OR
                (
                  U."VERSION"    IS NULL
                AND L.IN_VERSION IS NULL
                )
              )
            )
          THEN 'N'
          ELSE 'Y'
        END "DELETED",
        CASE
          WHEN
            (
              U.UPI = L.COMPARABLE_PROT_UPI
            AND
              (
                U."VERSION" = L.IN_VERSION
              OR
                (
                  U."VERSION"    IS NULL
                AND L.IN_VERSION IS NULL
                )
              )
            )
          THEN NVL (L.IN_TAXID, U.TAXID)
          ELSE U.TAXID
        END "TAXID"
      FROM
        load_retro_tmp l
      WHERE
        u.deleted                = 'N'
      AND u.AC                   = l.IN_AC
      AND U.DBID                 = L.IN_DBID
      AND U.last                 < L.IN_LOAD_RELEASE
      AND l.COMPARABLE_PROT_UPI IS NOT NULL
    )
  WHERE
    U.DELETED = 'N'
  AND EXISTS
    (
      SELECT
        /*+ PARALLEL (X 4) */
        1
      FROM
        LOAD_RETRO_TMP x
      WHERE
        u.AC                     = x.IN_AC
      AND U.DBID                 = X.IN_DBID
      AND U.last                 < X.IN_LOAD_RELEASE
      AND x.COMPARABLE_PROT_UPI IS NOT NULL
    );

  -- inserts new xref
  -- for each deleted = 'Y'
  -- in the above update
  INSERT
    /*+ PARALLEL (XREF 4) */
  INTO
    RNACEN.xref
    (
      ac,
      dbid,
      version,
      version_i,
      upi,
      created,
      last,
      deleted,
      taxid
    )
  SELECT
    /*+ PARALLEL (L 4) PARALLEL (X 4) */
    l.IN_AC,
    l.IN_DBID,
    l.IN_VERSION,
    CASE
      WHEN
        (
          X.UPI != L.COMPARABLE_PROT_UPI
        )
      OR
        (
          X."VERSION" != L.IN_VERSION
        )
      OR
        (
          X."VERSION"    IS NULL
        AND L.IN_VERSION IS NOT NULL
        )
      OR
        (
          X."VERSION"    IS NOT NULL
        AND L.IN_VERSION IS NULL
        )
      THEN X.VERSION_I + 1
      ELSE X.VERSION_I
    END VERSION_I_NEW,
    COMPARABLE_PROT_UPI,
    IN_LOAD_RELEASE CREATED,
    IN_LOAD_RELEASE LAST,
    'N' DELETED,
    IN_TAXID
  FROM
    RNACEN.XREF X,
    LOAD_RETRO_TMP L
  WHERE
    X.AC        = L.IN_AC
  AND X.DBID    = L.IN_DBID
  AND X.DELETED = 'Y'
  AND x.last    = v_previous_release
    -- instead of
    -- AND x.last                 < l.IN_LOAD_RELEASE
  AND l.COMPARABLE_PROT_UPI IS NOT NULL
  AND NOT
    (
      X.UPI = L.COMPARABLE_PROT_UPI
    AND
      (
        X."VERSION" = L.IN_VERSION
      OR
        (
          X."VERSION"    IS NULL
        AND L.IN_VERSION IS NULL
        )
      )
    );

  -- inserts a new xref
  -- for dbid, ac where xref:
  -- a) does not exists at all (VERSION_I = 1)
  -- b) does not exists because has been deleted
  --    1) and was having the same UPI, Version (VERSION_I = OLD_VERSION_I)
  --    2) and was having a different UPI, VERSION (VERSION_I = OLD_VERSION_I
  -- + 1)
  INSERT
    /*+ PARALLEL (XREF 4) */
  INTO
    RNACEN.xref
    (
      ac,
      dbid,
      version,
      version_i,
      upi,
      created,
      last,
      deleted,
      taxid
    )
  SELECT
    /*+ PARALLEL (T 4) */
    T.IN_AC,
    T.IN_DBID,
    T.IN_VERSION,
    CASE T.COMPARABLE_PROT_UPI
      WHEN
        (
          SELECT
            /*+ PARALLEL (PREVIOUS_XREF 4) */
            UNIQUE(PREVIOUS_XREF.UPI)
          FROM
            RNACEN.XREF PREVIOUS_XREF
          WHERE
            PREVIOUS_XREF.AC          = T.IN_AC
          AND PREVIOUS_XREF.VERSION_I = T.MAX_PREVIOUS_XREF_VERSION_I
          AND PREVIOUS_XREF.DBID      = T.IN_DBID
        )
      THEN T.MAX_PREVIOUS_XREF_VERSION_I
      ELSE T.MAX_PREVIOUS_XREF_VERSION_I + 1
    END VERSION_I,
    T.COMPARABLE_PROT_UPI UPI,
    T.IN_LOAD_RELEASE CREATED,
    T.IN_LOAD_RELEASE LAST,
    'N' DELETED,
    T.IN_TAXID TAXID
  FROM
    (
      SELECT
        /*+ PARALLEL (L 4) */
        L.IN_AC,
        L.IN_DBID,
        L.IN_VERSION,
        NVL (
        (
          SELECT
            /*+ PARALLEL (PREVIOUS_XREF 4) */
            MAX(PREVIOUS_XREF.VERSION_I)
          FROM
            RNACEN.XREF PREVIOUS_XREF
          WHERE
            PREVIOUS_XREF.AC     = L.IN_AC
          AND PREVIOUS_XREF.DBID = L.IN_DBID
        )
        , 0 ) MAX_PREVIOUS_XREF_VERSION_I,
        L.COMPARABLE_PROT_UPI,
        L.IN_LOAD_RELEASE,
        L.IN_TAXID
      FROM
        LOAD_RETRO_TMP L
      WHERE
        L.COMPARABLE_PROT_UPI IS NOT NULL
      AND NOT EXISTS
        (
          SELECT
            /*+ PARALLEL (X 4) */
            1
          FROM
            RNACEN.XREF X
          WHERE
            X.AC        = L.IN_AC
          AND X.DBID    = L.IN_DBID
          AND X.DELETED = 'N'
          AND X.LAST   <= L.IN_LOAD_RELEASE
            -- using also = in X.LAST <= L.IN_LOAD_RELEASE
            -- to eliminate from the selection
            -- the xref inserted right above after updating
            -- if active xref exists with another protein/version
            -- MUST have been updated/inserted above
            -- nothing else to insert here
        )
    )
    T;

  COMMIT;

  EXECUTE IMMEDIATE 'ALTER SESSION FORCE PARALLEL DML';

END incremental_release;

/*
* For a given database release, the rows of its corresponding 'load_X' table
* are loaded.
*/
FUNCTION PROCESS_LOAD_TABLE_SETBASED(
    p_in_dbid         IN RNACEN.RNC_DATABASE.ID%TYPE,
    p_in_load_release IN RNACEN.RNC_RELEASE.ID%TYPE )
  RETURN PLS_INTEGER
IS

  -- Anton: the following two declarations don't seem to be used
  TYPE t_load_record
  IS
    RECORD
    (
      crc64 RNACEN.rna.crc64%TYPE,
      LEN RNACEN.rna.len%TYPE,
      seq_short RNACEN.rna.seq_short%TYPE,
      seq_long RNACEN.rna.seq_long%TYPE,
      ac RNACEN.xref.ac%TYPE,
      version RNACEN.xref.version%TYPE,
      md5 RNACEN.rna.md5%TYPE);

  TYPE t_load_cursor
  IS
    REF
    CURSOR;
      v_load t_load_record;
      c_load t_load_cursor;
      v_previous_release RNACEN.rnc_release.id%TYPE;
      v_release_type RNACEN.rnc_release.release_type%TYPE;
      v_previous_release_type RNACEN.rnc_release.release_type%TYPE;
      v_upi RNACEN.xref.upi%TYPE;

  l_count PLS_INTEGER := 0;

BEGIN

  DBMS_OUTPUT.PUT_LINE('Loading release: ' || p_in_load_release);

  v_release_type     := RNACEN.release.get_release_type(p_in_load_release);
  v_previous_release := RNACEN.release.get_previous_release(p_in_dbid, p_in_load_release);

  IF v_previous_release     IS NOT NULL THEN
    v_previous_release_type := RNACEN.release.get_release_type(v_previous_release);
  END IF;

  -------------------------
  -- loading temp tables --
  -------------------------
  load_retro_tmp_table(p_in_dbid, p_in_load_release);
  load_md5_stats_table;
  load_md5_collisions_table;
  load_md5_new_sequences_table;
  set_comparable_prot_upi;
  -- store new sequences, no work on xrefs yet
  store_new_sequences;

  -- begin xref analysis
  IF v_release_type = 'F'
  THEN

    prepare_xref_pel_tables;
    populate_pel_tables1(v_previous_release);
    populate_pel_tables2;

    --code_to_split(p_in_dbid, v_previous_release, v_release_type);
    load_upi_max_versions_table(p_in_dbid);
    load_max_versions_table;
    populate_pel_tables3(p_in_dbid);

    /*
    * For full releases, retire those entries which doesn't appear anymore
    * in the latest release.
    *
    * Partition Exchange Loading (PEL) is better suited for
    * (F)ull and (C)umulative releases, where the current data replace totally
    * the previous ones
    * but incremental logic will be managed also, just in case
    *
    */
    IF v_release_type = 'F' AND v_previous_release IS NOT NULL THEN
      populate_pel_tables4(p_in_dbid, v_previous_release);
    -- Anton: this condition never seems to be true if only F releases are considered
    ELSIF V_PREVIOUS_RELEASE IS NOT NULL THEN
      populate_pel_tables5(p_in_dbid, v_previous_release, v_release_type);
    END IF;

    -- Anton: NB! don't forget to uncomment the next line
    prepare_pel_exchange;

    -- main PEL statements
     EXECUTE IMMEDIATE 'ALTER TABLE XREF EXCHANGE SUBPARTITION XREF_P' || TO_CHAR ( p_in_dbid ) ||
       '_NOT_DELETED WITH TABLE XREF_PEL_NOT_DELETED INCLUDING INDEXES WITHOUT VALIDATION';
     EXECUTE IMMEDIATE 'ALTER TABLE XREF EXCHANGE SUBPARTITION XREF_P' || TO_CHAR ( p_in_dbid ) ||
       '_DELETED WITH TABLE XREF_PEL_DELETED INCLUDING INDEXES WITHOUT VALIDATION';

  ELSE

    incremental_release(v_previous_release, p_in_dbid);

  END IF;

  /*
  * Set release status as 'Done' in RNC_release and update 'current_release' in
  * RNC_database.
  */
  EXECUTE IMMEDIATE 'ALTER SESSION DISABLE PARALLEL DML';
  -- the next statement causes ORA-12838 without disabling parallel DML
  RNACEN.release.set_release_status( p_in_load_release, 'D' );
  RNACEN.DATABASE.set_current_release( p_in_dbid, p_in_load_release );

  RETURN l_count;

END process_load_table_setbased;

/*
* Performs the load of a given database. It first processes all rows
* in LOAD_XXX table, then it checks that the load didn't delete all
* records from the previous release, and then it sends an email notifing the
load.
*/
PROCEDURE load_database_release_setbased(
    in_dbid         IN RNACEN.rnc_database.id%TYPE,
    in_load_release IN RNACEN.rnc_release.id%TYPE)
IS
  v_count PLS_INTEGER;
BEGIN

  v_count := process_load_table_setbased( in_dbid, in_load_release );

  -- restore logging
  -- two extra pel exchange statements here

END load_database_release_setbased;

/*
* This procedure is called by a scheduled Oracle job and runs every 5 minutes.
* This procedure executes 'load_database_release_setbased' for each release
* with status 'L' i.e. a release awaiting to be loaded.
*/
PROCEDURE load_job_setbased
IS

  CURSOR c_load
  IS
    SELECT
      dbid,
      id,
      release_type,
      release_date,
      force_load
    FROM
      RNACEN.rnc_release
    WHERE
      status = 'L'
    ORDER BY
      id;

BEGIN

  DBMS_OUTPUT.PUT_LINE('load_job_setbased');

  FOR v_load IN c_load
  LOOP

    DBMS_OUTPUT.PUT_LINE('Release type = ' || v_load.release_type);

    load_database_release_setbased (
      in_dbid         => v_load.dbid,
      in_load_release => v_load.id
    );

  END LOOP;

END load_job_setbased;

END load_retro_setbased;