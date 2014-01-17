-- Copyright [2009-2014] EMBL-European Bioinformatics Institute
--
-- Licensed under the Apache License, Version 2.0 (the "License");
-- you may not use this file except in compliance with the License.
-- You may obtain a copy of the License at
--
--      http://www.apache.org/licenses/LICENSE-2.0
--
-- Unless required by applicable law or agreed to in writing, software
-- distributed under the License is distributed on an "AS IS" BASIS,
-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
-- See the License for the specific language governing permissions and
-- limitations under the License.

set define off

create or replace package expunge is

  /*
   * Expunges the first or last release for the given database.
   */
   procedure expunge( in_release in RNACEN.rnc_release.id%TYPE );

   pragma restrict_references (expunge, RNPS, WNPS);

end;
/
create or replace package body expunge is




   procedure expunge(
      in_release in RNACEN.rnc_release.id%TYPE )
   is
     v_dbid RNACEN.rnc_database.id%TYPE;

     v_previous_release RNACEN.rnc_release.id%TYPE;
     v_next_release RNACEN.rnc_release.id%TYPE;

   begin

     select dbid into v_dbid

        from rnc_release where id = in_release;


     if v_dbid = null then
        return;
     end if;

     v_previous_release := RNACEN.release.get_previous_release( v_dbid, in_release );
     v_next_release     := RNACEN.release.get_next_release    ( v_dbid, in_release );

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

           RNACEN.database.set_current_release(v_dbid, NULL);


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

           RNACEN.database.set_current_release(v_dbid, v_previous_release);

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
end;
/
set define on
