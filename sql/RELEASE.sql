set define off

create or replace package release
is
   /*
    * The release does not exists.
    */
   release_not_found exception;

   /*
    * Returns the current release given the database identifier. If
    * the database identifier does not exists then returns null.
    */
   function get_current_release(
      in_dbid in RNACEN.RNC_database.id%TYPE)

         return RNACEN.rnc_release.id%TYPE;

   /*
    * Returns the load release given the database identifier. If
    * the database identifier does not exists then returns null.
    */
   function get_load_release(
      in_dbid in RNACEN.RNC_database.id%TYPE)
         return RNACEN.rnc_release.id%TYPE;

   /*
    * Returns the previous release given the database identifier
    * and release identifier. If the database or release identifier

    * does not exists then returns null.
    */
   function get_previous_release(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_id in RNACEN.rnc_release.id%TYPE)
         return RNACEN.rnc_release.id%TYPE;


   /*
    * Returns the next release given the database identifier
    * and release identifier. If the database or release identifier
    * does not exists then returns null.
    */

   function get_next_release(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_id in RNACEN.rnc_release.id%TYPE)
         return RNACEN.rnc_release.id%TYPE;

   /*
    * Returns the latest release given the database identifier. If
    * the database identifier does not exists then returns null.
    */
   function get_latest_release(
      in_dbid in RNACEN.RNC_database.id%TYPE)
         return RNACEN.rnc_release.id%TYPE;


   /*
    * Returns the release type given the release identifier. If
    * the release identifier does not exists then returns null.
    */
   function get_release_type(
      in_release_id in RNACEN.rnc_release.id%TYPE)
         return RNACEN.rnc_release.release_type%TYPE;

   /*
    * Returns the release status given the release identifier. If
    * the release identifier does not exists then returns null.
    */
   function get_release_status(

      in_release_id in RNACEN.rnc_release.id%TYPE)
         return RNACEN.rnc_release.status%TYPE;
   /*
    * Returns the release identifier given the database identifier and
    * the release date. If the database identifier or the release
    * date does not exists then returns null.
    */
   function get_release_id(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_date in RNACEN.rnc_release.release_date%TYPE)
         return RNACEN.rnc_release.id%TYPE;

   /*

    * Sets the release status given the release identifier and
    * the release status.
    */
   procedure set_release_status(
      in_release_id in RNACEN.rnc_release.id%TYPE,
      release_status in RNACEN.rnc_release.status%TYPE);

   /*
    * Returns the checkpoint, which is the latest imported release
    */
   function get_checkpoint
      return RNACEN.rnc_release.id%TYPE;


   /*
    * Given a release, this function returns the number of entries
    * that were retired when loading that release.
    */
   function get_retired_count(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release in RNACEN.rnc_release.id%TYPE )
   return PLS_INTEGER;

   /*
    * Given a release, this function returns the number of active entries
    * of that release.
    */

   function get_active_count(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release in RNACEN.rnc_release.id%TYPE )
   return PLS_INTEGER;


   pragma restrict_references (get_current_release, WNDS, RNPS, WNPS);
   pragma restrict_references (get_load_release, WNDS, RNPS, WNPS);
   pragma restrict_references (get_previous_release, WNDS, RNPS, WNPS);
   pragma restrict_references (get_next_release, WNDS, RNPS, WNPS);
   pragma restrict_references (get_latest_release, WNDS, RNPS, WNPS);
   pragma restrict_references (get_release_type, WNDS, RNPS, WNPS);
   pragma restrict_references (get_release_status, WNDS, RNPS, WNPS);

   pragma restrict_references (get_release_id, WNDS, RNPS, WNPS);
   pragma restrict_references (set_release_status, RNPS, WNPS);
   pragma restrict_references (get_checkpoint, RNPS, WNPS);

end release;
/
create or replace package body release

is

   function get_current_release(
      in_dbid in RNACEN.RNC_database.id%TYPE)
         return RNACEN.rnc_release.id%TYPE

   is
      v_id RNACEN.rnc_release.id%TYPE;
   begin
      select current_release into v_id from RNACEN.RNC_database
         where id = in_dbid;

      return v_id;
   exception
      when no_data_found then
         return v_id;
   end;

   function get_load_release(

      in_dbid in RNACEN.RNC_database.id%TYPE)
         return RNACEN.rnc_release.id%TYPE
   is
      v_id RNACEN.rnc_release.id%TYPE;
   begin
      select id into v_id from RNACEN.rnc_release
         where dbid = in_dbid and status <> 'D';

      return v_id;
   exception
      when no_data_found then
         return v_id;
   end;


   function get_previous_release(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_id in RNACEN.rnc_release.id%TYPE)
        return RNACEN.rnc_release.id%TYPE
   is
      v_id RNACEN.rnc_release.id%TYPE;
   begin
      select max (id) into v_id from RNACEN.rnc_release
         where dbid = in_dbid AND id < in_release_id;

      return v_id;
   exception

      when no_data_found then
         return v_id;
   end;

   function get_next_release(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_id in RNACEN.rnc_release.id%TYPE)
        return RNACEN.rnc_release.id%TYPE
   is
      v_id RNACEN.rnc_release.id%TYPE;
   begin
      select min (id) into v_id from RNACEN.rnc_release
         where dbid = in_dbid AND id > in_release_id;


      return v_id;
   exception
      when no_data_found then
         return v_id;
   end;

   function get_latest_release(
      in_dbid in RNACEN.RNC_database.id%TYPE)
         return RNACEN.rnc_release.id%TYPE
   is
      v_id RNACEN.rnc_release.id%TYPE;
   begin

      select max (id) into v_id from RNACEN.rnc_release
         where dbid = in_dbid;

      return v_id;
   exception
      when no_data_found then
         return v_id;
   end;

   function get_release_type(
      in_release_id in RNACEN.rnc_release.id%TYPE)
         return RNACEN.rnc_release.release_type%TYPE
   is

      v_release_type RNACEN.rnc_release.release_type%TYPE;
   begin
      select release_type into v_release_type
         from RNACEN.rnc_release where id = in_release_id;

      return v_release_type;
   exception
      when no_data_found then
         raise release_not_found;
   end;

   function get_release_status(
      in_release_id in RNACEN.rnc_release.id%TYPE)

         return RNACEN.rnc_release.status%TYPE
   is
      v_release_status RNACEN.rnc_release.status%TYPE;
   begin
      select status into v_release_status
         from RNACEN.rnc_release where id = in_release_id;

      return v_release_status;
   exception
      when no_data_found then
         raise release_not_found;
   end;


   function get_release_id(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_date in RNACEN.rnc_release.release_date%TYPE)
         return RNACEN.rnc_release.id%TYPE
   is
      v_release_id RNACEN.rnc_release.id%TYPE;
   begin
      select id into v_release_id
         from RNACEN.rnc_release where dbid = in_dbid and
            release_date = in_release_date;

      return v_release_id;
   exception

      when no_data_found then
         return v_release_id;
   end;

   procedure set_release_status(
      in_release_id in RNACEN.rnc_release.id%TYPE,
      release_status in RNACEN.rnc_release.status%TYPE
      )
   is
   begin
      update RNACEN.rnc_release set status = release_status
         where id = in_release_id;
   end;


  function get_checkpoint
     return RNACEN.rnc_release.id%TYPE
  is
     v_checkpoint RNACEN.rnc_release.id%TYPE;
  begin
     select max(id) into v_checkpoint
        from RNACEN.rnc_release where status = 'D';

     return v_checkpoint;
  end;



/*
 * get_retired_count
 */

function get_retired_count(
   in_dbid in RNACEN.RNC_database.id%TYPE,
   in_release in RNACEN.rnc_release.id%TYPE )
   return PLS_INTEGER
is

   v_count PLS_INTEGER := 0;

   v_previous RNACEN.rnc_release.id%TYPE;


begin

   -- get id of previous database's release

    v_previous := get_previous_release(in_dbid, in_release);

   -- get the number of entries retired during the in_release

   select count(*) into v_count
      from RNACEN.xref
      where dbid = in_dbid and
            created < in_release and -- logically redundant but it optimizes the query

            last = v_previous;

   return v_count;

end;


/*
 * get_active_count
 */

function get_active_count(
   in_dbid in RNACEN.RNC_database.id%TYPE,

   in_release in RNACEN.rnc_release.id%TYPE )
   return PLS_INTEGER
is

   v_count PLS_INTEGER := 0;

begin

   -- get the number of active entries in the in_release

   select count(*) into v_count
      from RNACEN.xref
      where dbid = in_dbid and

            created <= in_release and
            ( deleted = 'N' or
              last >= in_release );

   return v_count;

end;


end release;
/
set define on
