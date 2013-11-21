set define off

create or replace function         get_upi_table(
   in_string in varchar2)
      return type_upi_table
   as
      v_string long := in_string || ',';

      v_pos number;

      v_upi_table type_upi_table := RNACEN.type_upi_table();

      begin
         loop
            v_pos := instr ( v_string, ',' );


            exit when ( nvl ( v_pos, 0 ) = 0 );

            v_upi_table.extend;

            v_upi_table ( v_upi_table.count ) := trim (
               substr ( v_string, 1, v_pos - 1 ) );

            v_string := substr( v_string, v_pos + 1);

         end loop;
         return (v_upi_table);
    end;
/
set define on
