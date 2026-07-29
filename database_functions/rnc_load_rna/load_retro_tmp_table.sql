CREATE OR REPLACE FUNCTION rnc_load_rna.load_retro_tmp_table(p_in_dbid bigint, p_in_load_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    EXECUTE 'drop index if exists load_retro_tmp$ac_dbid_upi';

    EXECUTE 'truncate table load_retro_tmp';

/*
    EXECUTE 'INSERT INTO LOAD_RETRO_TMP (
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
      $1 p_in_dbid,
      $2 p_in_load_release,
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
              AND in_seq_long = prot_seq_long
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
      END
    FROM
      (
      WITH
        -- distinct_loaded_rows contains sequences from the staging table
        distinct_loaded_rows AS (
          SELECT
            DISTINCT CRC64,
            LEN,
            SEQ_SHORT,
            seq_long MY_SEQ_LONG,
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
        l.MY_SEQ_LONG AS in_seq_long,
        l.ac in_ac,
        l.version in_version,
        L.MD5 IN_MD5,
        l.taxid in_taxid,
        p.id prot_id,
        p.seq_short prot_seq_short,
        p.seq_long prot_seq_long,
        p.len prot_len,
        P.UPI PROT_UPI
      FROM rnacen.rna p RIGHT OUTER JOIN distinct_loaded_rows l ON (p.md5 = l.md5) ) alias1' USING p_in_dbid, p_in_load_release;
*/

    EXECUTE 'drop table if exists distinct_loaded_rows_tmp';

    EXECUTE '
      create temporary table distinct_loaded_rows_tmp 
          on commit preserve rows 
      as
      select distinct crc64, len, seq_short, seq_long my_seq_long, ac, version, taxid, md5
        from load_rnacentral
    ';
                        
    EXECUTE 'analyze distinct_loaded_rows_tmp';

    insert into load_retro_tmp (
        in_dbid ,
        in_load_release ,
        in_crc64 ,
        in_len ,
        in_seq_short ,
        in_seq_long ,
        in_ac ,
        in_version ,
        in_md5 ,
        in_taxid ,
        comparable_prot_upi
      )
    select
      p_in_dbid, 
      p_in_load_release,
      in_crc64,
      in_len,
      in_seq_short,
      in_seq_long,
      in_ac,
      in_version,
      in_md5,
      in_taxid,
      case
        when prot_upi is not null
        then -- a rna with the same md5 exists
          case
            when
              (
                (
                  in_len   > 4000
                and in_len = prot_len
                )
              and in_seq_long = prot_seq_long
              )
            then -- a rna with the same md5, len, sequence (lob) exists
              prot_upi
            when
              (
                (
                  in_len  <= 4000
                and in_len = prot_len
                )
              and in_seq_short = prot_seq_short
              )
            then -- a rna with the same md5, len, sequence (not lob) exists
              prot_upi
            else -- the rna is not really the same
              null
          end
      end
    from
      (
      select
        l.crc64 in_crc64,
        l.len in_len,
        l.seq_short in_seq_short,
        l.my_seq_long as in_seq_long,
        l.ac in_ac,
        l.version in_version,
        l.md5 in_md5,
        l.taxid in_taxid,
        p.id prot_id,
        p.seq_short prot_seq_short,
        p.seq_long prot_seq_long,
        p.len prot_len,
        p.upi prot_upi
      from rna p right outer join distinct_loaded_rows_tmp l on (p.md5 = l.md5) ) alias1;

    --COMMIT;

    execute 'analyze load_retro_tmp';

  END;

$function$

