set define off

create or replace public class Upi
{
  public static String getUpi( long id )
  {
     String str = Long.toHexString( id ).toUpperCase();
     return "UPI0000000000".substring(0, 13 - str.length() ) + str;
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
