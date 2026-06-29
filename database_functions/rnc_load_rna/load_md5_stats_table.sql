CREATE OR REPLACE FUNCTION rnc_load_rna.load_md5_stats_table()
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
BEGIN

    EXECUTE 'TRUNCATE TABLE load_md5_stats';
    INSERT INTO RNACEN.LOAD_MD5_STATS (
        IN_MD5,
        CNT,
        CNT_DST_SEQ_SHORT,
        CNT_DST_SEQ_LONG
      )
    SELECT
      IN_MD5,
      COUNT(*) cnt,
      COUNT(DISTINCT IN_SEQ_SHORT) CNT_DST_SEQ_SHORT,
      COUNT(DISTINCT IN_SEQ_LONG) CNT_DST_SEQ_LONG
    FROM (
        SELECT
          l.in_seq_short,
          l.in_seq_long,
          l.in_md5
        FROM
          LOAD_RETRO_TMP L
        WHERE l.comparable_prot_upi IS NULL
      ) alias5
    GROUP BY
      IN_MD5;

    --COMMIT;

  END;

$function$

