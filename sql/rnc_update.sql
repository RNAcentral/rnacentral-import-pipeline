create or replace
PACKAGE BODY RNC_UPDATE AS

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
    v_previous_release RNACEN.rnc_release.id%TYPE;
  BEGIN

    DBMS_OUTPUT.PUT_LINE('Loading release: ' || p_in_load_release);

    -- initial logging
    RNC_LOGGING.log_release_start(p_in_dbid, p_in_load_release);

    -- load sequences
    RNC_LOAD_RNA.load_rna(p_in_dbid, p_in_load_release);

    -- load xrefs
    v_previous_release := RNACEN.release.get_previous_release(p_in_dbid, p_in_load_release);
    RNC_LOAD_XREF.load_xref(v_previous_release, p_in_dbid);

    mark_as_done(p_in_dbid, p_in_load_release);

    RNC_LOGGING.log_release_end(p_in_dbid, p_in_load_release, v_previous_release);

  END load_release;

  /*
    Iterates over all releases with status 'L' and load them into the database.
  */
  PROCEDURE new_update
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

    DBMS_OUTPUT.put_line('Launching an update');

    FOR v_load IN c_load
    LOOP
      load_release(
        p_in_dbid         => v_load.dbid,
        P_in_load_release => v_load.id
      );
    END LOOP;

  END new_update;

END RNC_UPDATE;