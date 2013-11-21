set define off

create or replace PACKAGE RNC_UPDATE AS

  /*
    This is the main entry point for importing new data into the RNAcentral database.
  */
  PROCEDURE new_update(p_release_type IN RNACEN.rnc_release.release_type%TYPE DEFAULT 'F');

  PROCEDURE update_literature_references;

  PROCEDURE update_accession_info;

  PROCEDURE update_composite_ids;


  PROCEDURE update_rnc_accessions;

END RNC_UPDATE;
/
create or replace PACKAGE BODY RNC_UPDATE AS


  /*
  * Based on the rnc_ac_info and rnc_composite_ids
  * update a denormalized table combining the other two tables.
  * Useful for convenient access using Django ORM.
  */
  PROCEDURE update_rnc_accessions

  IS
  BEGIN

    DBMS_OUTPUT.put_line('Updating rnc_accessions');

    -- merge regular ids from rnc_ac_info
    MERGE /*+ parallel */ INTO rnc_accessions t1
    USING (SELECT * FROM rnc_ac_info) t2
    ON (t1.accession = t2.ac)
    WHEN MATCHED THEN UPDATE SET
        t1.parent_ac = t2.parent_ac,
        t1.seq_version = t2.seq_version,
        t1.feature_start = t2.feature_start,

        t1.feature_end = t2.feature_end,
        t1.feature_name = t2.feature_name,
        t1.ordinal = t2.ordinal,
        t1.division = t2.division,
        t1.keywords = t2.keywords,
        t1.description = t2.description,
        t1.species = t2.species,
        t1.organelle = t2.organelle,
        t1.classification = t2.classification,
        t1.project = t2.project
    WHEN NOT MATCHED THEN INSERT
        (t1.accession,
         t1.parent_ac,

         t1.seq_version,
         t1.feature_start,
         t1.feature_end,
         t1.feature_name,
         t1.ordinal,
         t1.division,
         t1.keywords,
         t1.description,
         t1.species,
         t1.organelle,
         t1.classification,
         t1.project,
         t1.is_composite)

    VALUES
        (t2.ac,
         t2.parent_ac,
         t2.seq_version,
         t2.feature_start,
         t2.feature_end,
         t2.feature_name,
         t2.ordinal,
         t2.division,
         t2.keywords,
         t2.description,
         t2.species,
         t2.organelle,

         t2.classification,
         t2.project,
         'N');

    -- merge composite ids from rnc_composite_ids
    MERGE /*+ parallel */ INTO rnc_accessions t1
    USING (SELECT distinct t3.*, t4.COMPOSITE_ID, t4.DATABASE, t4.OPTIONAL_ID, t4.EXTERNAL_ID FROM rnc_ac_info t3, rnc_composite_ids t4 WHERE t3.ac = t4.ac) t2
    ON (t1.accession = t2.COMPOSITE_ID)
    WHEN MATCHED THEN UPDATE SET
        t1.parent_ac = t2.parent_ac,
        t1.seq_version = t2.seq_version,
        t1.feature_start = t2.feature_start,
        t1.feature_end = t2.feature_end,

        t1.feature_name = t2.feature_name,
        t1.ordinal = t2.ordinal,
        t1.division = t2.division,
        t1.keywords = t2.keywords,
        t1.description = t2.description,
        t1.species = t2.species,
        t1.organelle = t2.organelle,
        t1.classification = t2.classification,
        t1.project = t2.project,
        t1.non_coding_id = t2.composite_id,
        t1.database = t2.database,
        t1.optional_id = t2.optional_id,
        t1.external_id = t2.external_id

    WHEN NOT MATCHED THEN INSERT
        (t1.accession,
         t1.parent_ac,
         t1.seq_version,
         t1.feature_start,
         t1.feature_end,
         t1.feature_name,
         t1.ordinal,
         t1.division,
         t1.keywords,
         t1.description,
         t1.species,
         t1.organelle,

         t1.classification,
         t1.project,
         t1.is_composite,
         t1.non_coding_id,
         t1.database,
         t1.external_id,
         t1.optional_id)
    VALUES
        (t2.composite_id,
         t2.parent_ac,
         t2.seq_version,
         t2.feature_start,
         t2.feature_end,

         t2.feature_name,
         t2.ordinal,
         t2.division,
         t2.keywords,
         t2.description,
         t2.species,
         t2.organelle,
         t2.classification,
         t2.project,
         'Y',
         t2.ac,
         t2.database,
         t2.external_id,

         t2.optional_id);

    COMMIT;

  END update_rnc_accessions;

  /*
  *
  */
  PROCEDURE update_composite_ids
  IS
  BEGIN


    DBMS_OUTPUT.put_line('Updating composite ids');

    MERGE INTO rnc_composite_ids t1
    --ignore duplicates
    USING (SELECT * FROM rnc_composite_ids_all WHERE ROWID IN (SELECT MIN (ROWID) FROM rnc_composite_ids_all GROUP BY ac)) t2
		ON (t1.composite_id = t2.composite_id)
		WHEN MATCHED THEN UPDATE SET
  		t1.OPTIONAL_ID = t2.OPTIONAL_ID
		WHEN NOT MATCHED THEN INSERT
		(
      t1.COMPOSITE_ID,
      t1.AC,
      t1.DATABASE,

      t1.OPTIONAL_ID,
      t1.EXTERNAL_ID
		)
		VALUES
		(
      t2.COMPOSITE_ID,
      t2.AC,
      t2.DATABASE,
      t2.OPTIONAL_ID,
      t2.EXTERNAL_ID
    );

    DBMS_OUTPUT.put_line('Composite ids updated');


  END update_composite_ids;


  /*
  * Merge in the accession-level information from the staging table.
  */
  PROCEDURE update_accession_info
  IS
  BEGIN

    DBMS_OUTPUT.put_line('Updating accession information');


    MERGE INTO rnc_ac_info t1
    --ignore duplicates
    USING (SELECT * FROM rnc_ac_info_all WHERE ROWID IN (SELECT MIN (ROWID) FROM rnc_ac_info_all GROUP BY ac)) t2
		ON (t1.ac = t2.ac)
		WHEN MATCHED THEN UPDATE SET
  		t1.DIVISION = t2.DIVISION,
    	t1.KEYWORDS = t2.KEYWORDS,
      t1.DESCRIPTION = t2.DESCRIPTION,
      t1.ORGANELLE = t2.ORGANELLE,
      t1.SPECIES = t2.SPECIES,
      t1.CLASSIFICATION = t2.CLASSIFICATION,
      t1."PROJECT" = t2."PROJECT"
		WHEN NOT MATCHED THEN INSERT

		(
  		t1.AC,
			t1.PARENT_AC,
			t1.SEQ_VERSION,
			t1.FEATURE_START,
			t1.FEATURE_END,
			t1.FEATURE_NAME,
			t1.ORDINAL,
			t1.DIVISION,
      t1."PROJECT",
			t1.KEYWORDS,
			t1.DESCRIPTION,
			t1.ORGANELLE,

			t1.SPECIES,
			t1.CLASSIFICATION
		)
		VALUES
		(
  		t2.AC,
			t2.PARENT_AC,
			t2.SEQ_VERSION,
			t2.FEATURE_START,
			t2.FEATURE_END,
			t2.FEATURE_NAME,
			t2.ORDINAL,
			t2.DIVISION,

      t2."PROJECT",
			t2.KEYWORDS,
			t2.DESCRIPTION,
			t2.ORGANELLE,
			t2.SPECIES,
			t2.CLASSIFICATION
		);

    DBMS_OUTPUT.put_line('Accession information updated');

  END update_accession_info;



  /*
  * Insert non-redundant literature references from the staging table into
  * the main table.
  */
  PROCEDURE update_literature_references
  IS
  BEGIN

    INSERT INTO rnc_references
    (
      md5,
      authors_md5,
      location,

      authors,
      title,
      pubmed,
      doi,
      publisher,
      editors
    )
    SELECT
      /*+ PARALLEL */
      in_md5,
      in_authors_md5,
      in_location,
      in_my_authors,

      in_title,
      in_pubmed,
      in_doi,
      in_publisher,
      in_editors
    FROM
      (
      WITH
        distinct_new_refs AS
        (
          SELECT DISTINCT
            md5,
            authors_md5,

            location,
            useful_clob_object(authors) MY_AUTHORS,
            title,
            pubmed,
            doi,
            publisher,
            editors
          FROM
            rnc_references_all
         )
      SELECT
        l.md5 in_md5,
        l.authors_md5 in_authors_md5,

        l.location in_location,
        TREAT(l.MY_AUTHORS AS USEFUL_CLOB_OBJECT).UCO AS in_my_authors,
        l.title in_title,
        l.pubmed in_pubmed,
        l.doi in_doi,
        l.publisher in_publisher,
        l.editors in_editors
      FROM
        distinct_new_refs l,
        rnc_references p
     WHERE p.md5 (+) = l.md5
        AND p.authors_md5 (+) = l.authors_md5
        and p.location (+) = l.location

        AND p.md5 IS NULL
        AND p.authors_md5 IS NULL
        and p.location is null
      );

    COMMIT;

  END update_literature_references;


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

    rnc_healthchecks.run_healthchecks();

  END new_update;

END RNC_UPDATE;
/
set define on
