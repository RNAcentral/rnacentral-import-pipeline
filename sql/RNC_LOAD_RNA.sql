-- Copyright [2009-2014] EMBL-European Bioinformatics Institute
-- 
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
-- 
--      http://www.apache.org/licenses/LICENSE-2.0
-- 
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

set define off

create or replace PACKAGE RNC_LOAD_RNA AS

  /*
    Package for populating the RNA table from the staging table LOAD_RNACENTRAL.
    Assign UPIs to all new sequences.
  */

  PROCEDURE load_rna( p_in_dbid IN RNACEN.rnc_database.ID%TYPE,
                      p_in_load_release IN RNACEN.rnc_release.ID%TYPE
                    );

END RNC_LOAD_RNA;
/

create or replace PACKAGE BODY RNC_LOAD_RNA AS


  /*
  * Loads the data from the staging table load_rnacentral
  * into a temporary table LOAD_RETRO_TMP.
  * Retrieves the UPI of the corresponding RNA
  * when the new and existing sequences are the same.
  */
  PROCEDURE load_retro_tmp_table(
      p_in_dbid         IN RNACEN.RNC_DATABASE.ID%TYPE,
      p_in_load_release IN RNACEN.rnc_release.ID%TYPE )
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
          SELECT
            DISTINCT CRC64,

            LEN,
            SEQ_SHORT,
            useful_clob_object(seq_long) MY_SEQ_LONG,
            AC,
            VERSION,
            TAXID,
            MD5
          FROM
            RNACEN.load_rnacentral
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
  * Maintaned for logging purposes only.
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
      CAST (RNACEN.UPI.GET_UPI(SEQ_UPI.CURRVAL) AS NVARCHAR2(13)) PROT_UPI

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

  * Deposit the new sequences into the main
  * table containing all sequences from all sources.
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

  /*
    The main procedure for importing data from the staging table
    into the main RNA table. Responsible for assigning UPIs.
  */

  PROCEDURE load_rna( p_in_dbid         IN RNACEN.rnc_database.ID%TYPE,
                      p_in_load_release IN RNACEN.rnc_release.ID%TYPE)
  AS
  BEGIN

    load_retro_tmp_table(p_in_dbid, p_in_load_release);
    load_md5_stats_table;
    load_md5_collisions_table;
    load_md5_new_sequences_table;
    set_comparable_prot_upi;
    store_new_sequences;

  END load_rna;


END RNC_LOAD_RNA;
/
set define on
