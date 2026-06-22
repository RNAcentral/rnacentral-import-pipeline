CREATE OR REPLACE FUNCTION rnc_update.prepare_releases(p_release_type character)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    q CURSOR
    FOR
      SELECT distinct
        d2.id
      FROM
          load_rnacentral_all d1,
          rnc_database d2
      WHERE
        d1.DATABASE = d2.descr;

    v_count_existing_releases numeric;

BEGIN

    SELECT count(*)
    INTO v_count_existing_releases
    FROM rnc_release
    WHERE status = 'L';

    IF (v_count_existing_releases > 0) THEN
      RAISE NOTICE 'Found releases to be loaded';
      RETURN;
    END IF;

    RAISE NOTICE 'Preparing the release table';

    FOR v_db IN q
    LOOP
      perform rnc_update.create_release(p_in_dbid => v_db.ID, p_release_type => p_release_type);
    END LOOP;

  END;

$function$

