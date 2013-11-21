set define off

create or replace package database
is

   /*
    * The database does not exists.
    */
   database_not_found exception;

   /*
    * Returns the database description given the database identifier. If
    * the given database identifier does not exists then returns null.
    */
   FUNCTION get_database_descr(

      in_dbid in RNACEN.RNC_database.id%TYPE)
         return RNACEN.RNC_database.descr%TYPE;

   /*
    * Returns the database identifier given the database description. If
    * the given database description does not exists then returns null.
    */
   function get_database_id(
      in_descr in RNACEN.RNC_database.descr%TYPE)
         return RNACEN.RNC_database.id%TYPE;

   /*
    * Sets the current release given the database identifier and the

    * load release identifier.
    */
   procedure set_current_release(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_id  in RNACEN.RNC_release.id%TYPE
      );

   pragma restrict_references (get_database_descr, WNDS, RNPS, WNPS);
   pragma restrict_references (get_database_id, WNDS, RNPS, WNPS);
   pragma restrict_references (set_current_release, RNPS, WNPS);

end database;
/

create or replace package body database

is
   function get_database_descr(
      in_dbid in RNACEN.RNC_database.id%TYPE)
         return RNACEN.RNC_database.descr%TYPE
   is
      v_descr RNACEN.RNC_database.descr%TYPE;
   begin
      select descr into v_descr from RNACEN.RNC_database
         where id = in_dbid;

      return v_descr;

   exception
      when no_data_found then
         return v_descr;
   end;

   function get_database_id(
      in_descr in RNACEN.RNC_database.descr%TYPE)
         return RNACEN.RNC_database.id%TYPE
   IS
      v_id RNACEN.RNC_database.id%TYPE;
   BEGIN
      select id into v_id from RNACEN.RNC_database
         where UPPER(descr) = UPPER(in_descr);


      return v_id;
   exception
      when no_data_found then
         return v_id;
   end;

   procedure set_current_release(
      in_dbid in RNACEN.RNC_database.id%TYPE,
      in_release_id in RNACEN.RNC_release.id%TYPE
      )
   is
   begin

      update RNC_database set current_release = in_release_id
         where id = in_dbid;
   end;

end database;
/
set define on
