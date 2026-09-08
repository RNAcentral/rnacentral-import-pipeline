CREATE OR REPLACE FUNCTION rnc_update.prepare_releases(p_release_type character)
 RETURNS void
 LANGUAGE plpgsql
 SECURITY DEFINER
AS $function$
DECLARE

    -- Bailing out whenever any release was pending let one abandoned load block
    -- every later one: no release got created for the staged database, and the
    -- stale pending release ran in its place.
    q CURSOR
    FOR
      SELECT distinct
        d2.id
      FROM
          load_rnacentral_all d1,
          rnc_database d2
      WHERE
        d1.DATABASE = d2.descr
      AND NOT EXISTS (
          SELECT
            1
          FROM
            rnc_release r
          WHERE
            r.dbid   = d2.id
          AND r.status = 'L'
        );

BEGIN

    RAISE NOTICE 'Preparing the release table';

    FOR v_db IN q
    LOOP
      perform rnc_update.create_release(p_in_dbid => v_db.ID, p_release_type => p_release_type);
    END LOOP;

  END;

$function$
