CREATE OR REPLACE FUNCTION expunge.expunge(in_release bigint)
 RETURNS void
 LANGUAGE plpgsql
 STABLE SECURITY DEFINER
AS $function$
DECLARE

     v_dbid RNACEN.rnc_database.id%TYPE;

     v_previous_release RNACEN.rnc_release.id%TYPE;
     v_next_release RNACEN.rnc_release.id%TYPE;


BEGIN

     select dbid into v_dbid

        from rnc_release where id = in_release;


     if v_dbid  IS NULL then
        return;
     end if;

     v_previous_release := release.get_previous_release( v_dbid, in_release );
     v_next_release     := release.get_next_release( v_dbid, in_release );

     if v_next_release     is     null and
        v_previous_release is     null then


          /*
           * Only release.
           */


          /* Created database cross-references.
           */
           delete from RNACEN.xref
           where
              dbid = v_dbid;

           perform database.set_current_release(v_dbid, NULL);


     elsif v_next_release  is     null and
        v_previous_release is not null then

          /*

           * Last release.
           */

          /* Created database cross-references.
           */
           delete from RNACEN.xref
           where
              dbid = v_dbid and

              created = in_release and
              last = in_release;

          /* Retired database cross-references.
           */

           update RNACEN.xref
              set deleted = 'N'
              where
                 dbid = v_dbid and
                 last = v_previous_release and
                 deleted = 'Y';


          /* Touched database cross-references.
           */
           update RNACEN.xref
              set last = v_previous_release
              where
                 dbid = v_dbid and

                 last = in_release;

           perform database.set_current_release(v_dbid, v_previous_release);

     elsif v_next_release  is not null and
        v_previous_release is     null then


          /*
           * First release.
           */

          /* Created database cross-references.
           */

           delete from RNACEN.xref
           where
              dbid = v_dbid and
              created = in_release and
              last = in_release;


          /* Touched database cross-references.
           */

          update RNACEN.xref
             set created = v_next_release
             where
                dbid = v_dbid and

                created = in_release;
     else
        return;
     end if;


     delete from xref_not_unique
        where dbid = v_dbid;

/*
 * The UNIRELEASE table has a FK to CV_RELEASE. It keeps track of the uniparc
 * release in the moment of the uniprot public release. If you are trying
 * to expunge a release used there, the procedure will fail.
 * You have first to manually amend the data there, after deciding what release

 * to mark there, instead of the one being deleted.
 */
     delete from rnc_release

        where id = in_release;
   end;

$function$
