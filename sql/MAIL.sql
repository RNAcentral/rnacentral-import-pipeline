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

create or replace package mail is

/*
 * Sends an email to one receipient.
 */
procedure mail (
   in_sender in varchar2,
   in_subject in varchar2,
   in_message in varchar2,
   in_receiver in varchar2
);

/*

 * Sends an email to two receipients.
 */
procedure mail (
   in_sender in varchar2,
   in_subject in varchar2,
   in_message in varchar2,
   in_receiver1 in varchar2,
   in_receiver2 in varchar2
);

end mail;
/
create or replace package body mail is



procedure mail (
   in_sender in varchar2,
   in_subject in varchar2,
   in_message in varchar2,
   in_receiver in varchar2
)
is
  v_serv varchar2(30)    := 'mailserv.ebi.ac.uk';
  v_crlf varchar2(2)     := CHR( 13 ) || CHR( 10 );
  v_port number          := 25;
  v_mesg varchar2(4000);

  v_conn utl_smtp.connection;
begin

  v_conn:= utl_smtp.open_connection( v_serv, v_port );

  utl_smtp.helo( v_conn, v_serv );
  utl_smtp.mail( v_conn, in_sender );
  utl_smtp.rcpt( v_conn, in_receiver );

  v_mesg :=
         'Date: '   || TO_CHAR( sysdate, 'dd Mon yy hh24:mi:ss' ) || v_crlf ||
         'From:'              || in_sender                        || v_crlf ||
         'Subject:'           || in_subject                       || v_crlf ||

         'To: '               || in_receiver                      || v_crlf ||
         ''         || v_crlf || in_message;


  utl_smtp.data( v_conn, v_mesg );
  utl_smtp.quit( v_conn );

end;

procedure mail (
   in_sender in varchar2,
   in_subject in varchar2,
   in_message in varchar2,

   in_receiver1 in varchar2,
   in_receiver2 in varchar2
)
is
  v_serv varchar2(30)    := 'mailserv.ebi.ac.uk';
  v_crlf varchar2(2)     := CHR( 13 ) || CHR( 10 );
  v_port number          := 25;
  v_mesg varchar2(4000);
  v_conn utl_smtp.connection;
begin

  v_conn:= utl_smtp.open_connection( v_serv, v_port );


  utl_smtp.helo( v_conn, v_serv );
  utl_smtp.mail( v_conn, in_sender );
  utl_smtp.rcpt( v_conn, in_receiver1 );
  utl_smtp.rcpt( v_conn, in_receiver2 );

  v_mesg :=
         'Date: '   || TO_CHAR( sysdate, 'dd Mon yy hh24:mi:ss' ) || v_crlf ||
         'From:'              || in_sender                        || v_crlf ||
         'Subject:'           || in_subject                       || v_crlf ||
         'To: '               || in_receiver1 || ';'
                              || in_receiver2                     || v_crlf ||
         ''         || v_crlf || in_message;



  utl_smtp.data( v_conn, v_mesg );
  utl_smtp.quit( v_conn );

end;

end mail;
/
set define on
