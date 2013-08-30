create or replace
PACKAGE BODY RNC_UPDATE AS

  /*
  * Move the data for the specified database into the staging table.
  */
  PROCEDURE move_staging_data (
    p_in_dbid IN RNACEN.rnc_database.ID%TYPE
  )
  IS
  BEGIN

    EXECUTE IMMEDIATE 'TRUNCATE TABLE load_rnacentral DROP STORAGE';

    INSERT INTO load_rnacentral (
      SELECT CRC64,
             LEN,
             SEQ_SHORT,
             SEQ_LONG,
             DATABASE,
             AC,
             OPTIONAL_ID,
             VERSION,
             TAXID,
             MD5
      FROM load_rnacentral_all d1, rnc_database d2
      WHERE d1.DATABASE = d2.descr AND d2.ID = p_in_dbid
    );

    COMMIT;

  END move_staging_data;

  /*
  * Create a new release for the specified database.
  */
  PROCEDURE create_release (
    p_in_dbid IN RNACEN.rnc_database.ID%TYPE,
    p_release_type IN RNACEN.rnc_release.release_type%TYPE
  )
  IS
    v_next_release INTEGER;
  BEGIN

    SELECT count(*) + 1 INTO v_next_release FROM rnc_release;

    DBMS_OUTPUT.put_line('Creating new release for database ' || TO_CHAR(p_in_dbid));

    INSERT INTO rnc_release
      (ID,
       dbid,
       release_date,
       release_type,
       status,
       rnc_release.TIMESTAMP,
       userstamp,
       descr,
       force_load)
    VALUES
      (
      v_next_release,
      p_in_dbid ,
      (SELECT to_char(trunc(SYSDATE),'dd-MON-yy') FROM dual),
      p_release_type,
      'L',
      (SELECT to_char(trunc(SYSDATE),'dd-MON-yy') FROM dual),
      'auto',
      '',
      'N'
    );

    COMMIT;

  END create_release;

  /*
  * Create new releases for all databases mentioned in the staging table.
  */
  PROCEDURE prepare_releases(
    p_release_type IN RNACEN.rnc_release.release_type%TYPE
  )
  IS
    CURSOR q
    IS
      SELECT distinct
        d2.id
      FROM
          load_rnacentral_all d1,
          rnc_database d2
      WHERE
        d1.DATABASE = d2.descr;
      v_count_existing_releases INTEGER;
  BEGIN

    SELECT count(*)
    INTO v_count_existing_releases
    FROM rnc_release
    WHERE status = 'L';

    IF (v_count_existing_releases > 0) THEN
      DBMS_OUTPUT.put_line('Found releases to be loaded');
      RETURN;
    END IF;

    DBMS_OUTPUT.put_line('Preparing the release table');

    FOR v_db IN q
    LOOP
      rnc_update.create_release(p_in_dbid => v_db.ID, p_release_type => p_release_type);
    END LOOP;

  END prepare_releases;

  /*
  * Set release status as 'Done' in rnc_release and update 'current_release' in
  * RNC_database.
  */
  PROCEDURE mark_as_done(
    p_in_dbid         IN RNACEN.rnc_database.ID%TYPE,
    p_in_load_release IN RNACEN.rnc_release.id%TYPE)
  IS
  BEGIN

    EXECUTE IMMEDIATE 'ALTER SESSION DISABLE PARALLEL DML';
    -- the next statement causes ORA-12838 without disabling parallel DML
    RNACEN.release.set_release_status( p_in_load_release, 'D' );
    RNACEN.DATABASE.set_current_release( p_in_dbid, p_in_load_release );
    COMMIT;
    EXECUTE IMMEDIATE 'ALTER SESSION ENABLE PARALLEL DML';

  END mark_as_done;

  /*
    Load one release into the database.
  */
  PROCEDURE load_release(
    p_in_dbid         IN RNACEN.rnc_database.id%TYPE,
    p_in_load_release IN RNACEN.rnc_release.id%TYPE)
  IS
    v_previous_release RNACEN.rnc_release.ID%TYPE;
    v_release_type RNACEN.rnc_release.release_type%TYPE;
  BEGIN

    DBMS_OUTPUT.PUT_LINE('Loading release: ' || p_in_load_release);

    -- initial logging
    RNC_LOGGING.log_release_start(p_in_dbid, p_in_load_release);

    -- load sequences
    RNC_LOAD_RNA.load_rna(p_in_dbid, p_in_load_release);

    -- load xrefs
    v_release_type := RNACEN.release.get_release_type(p_in_load_release);
    v_previous_release := RNACEN.release.get_previous_release(p_in_dbid, p_in_load_release);
    IF v_release_type = 'F' THEN
      RNC_LOAD_XREF.load_xref(v_previous_release, p_in_dbid);
    elsif v_release_type = 'I' THEN
      RNC_LOAD_XREF_INCREMENTAL.load_xref_incremental(v_previous_release, p_in_dbid);
    END IF;

    mark_as_done(p_in_dbid, p_in_load_release);

    RNC_LOGGING.log_release_end(p_in_dbid, p_in_load_release, v_previous_release);

    COMMIT;

  END load_release;

  /*
    Iterate over all releases with status 'L' and load them into the database.
  */
  PROCEDURE new_update (
    p_release_type IN RNACEN.rnc_release.release_type%TYPE DEFAULT 'F'
  )
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
        ID;
  BEGIN

    DBMS_OUTPUT.put_line('Launching an update');

    prepare_releases(p_release_type);

    FOR v_load IN c_load
    LOOP
      move_staging_data(p_in_dbid => v_load.dbid);
      load_release(
        p_in_dbid         => v_load.dbid,
        P_in_load_release => v_load.ID
      );
    END LOOP;

  END new_update;

END RNC_UPDATE;