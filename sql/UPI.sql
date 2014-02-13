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

create or replace public class Upi {

	public static String getUpi( long id ) {
		String str = Long.toHexString( id ).toUpperCase();
		return "URS0000000000".substring(0, 13 - str.length() ) + str;
	}

	public static void main (String args[]) {
		long l = Long.parseLong(args[0]);
		System.out.println(getUpi(l));
	}
}
create or replace package upi is


   /*
    * Returns a new protein identifier.

    */
   function get_upi (
      in_id in number)
         return varchar2;

   pragma restrict_references (get_upi, RNDS, WNDS, RNPS, WNPS);

end upi;
/

create or replace package body upi


is
   function get_upi (

      in_id in number)
         return varchar2
   as
      language java name 'Upi.getUpi(long) return java.lang.String';

end upi;
/
set define on
