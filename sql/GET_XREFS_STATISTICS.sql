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

create or replace PROCEDURE get_xrefs_statistics (
p_dbid number,
p_this_release number,
p_prev_release number := null
)
AS
l_dbid number := p_dbid;
l_this_release number := p_this_release;
L_PREV_RELEASE NUMBER := NVL (P_PREV_RELEASE, RNACEN.RELEASE.GET_PREVIOUS_RELEASE(L_DBID, L_THIS_RELEASE));
l_sql_stmt varchar2 (4000);
BEGIN

l_sql_stmt :=

'MERGE INTO RELEASE_STATS s
USING (
SELECT
:l_dbid dbid,
:l_this_release this_release,
:l_prev_release prev_release,
COUNT (*) ff_loaded_rows
FROM RNACEN.LOAD_' || RNACEN.DATABASE.get_database_descr( p_dbid ) || '
) q
on (s.this_release = q.this_release)
when matched then update
set
s.dbid = q.dbid,

s.prev_release = q.prev_release,
s.ff_loaded_rows = q.ff_loaded_rows
when not matched then insert
(
dbid,
this_release,
prev_release,
ff_loaded_rows)
values
(
q.dbid,
q.this_release,
q.prev_release,

q.ff_loaded_rows)';

-- DBMS_OUTPUT.PUT_LINE (L_SQL_STMT);
EXECUTE IMMEDIATE l_sql_stmt USING l_dbid, l_this_release, l_prev_release;

merge into RELEASE_STATS s
using (
select
l_dbid dbid,
l_this_release this_release,
l_prev_release prev_release,
retired_prev_releases,
retired_this_release,

retired_next_releases,
retired_total,
created_w_predecessors_v_1,
created_w_predecessors_v_gt1,
created_w_predecessors,
created_wo_predecessors_v_1,
created_wo_predecessors_v_gt1,
created_wo_predecessors,
active_created_prev_releases,
active_created_this_release,
active_created_next_releases,
created_this_release,
active_updated_this_release,

active_untouched_this_release,
active_total
from (
select
sum (retired_prev_releases) retired_prev_releases,
sum (retired_this_release) retired_this_release,
sum (retired_next_releases) retired_next_releases,
sum (retired_total) retired_total,
sum (case when version_i = 1 then created_w_predecessors else 0 end) created_w_predecessors_v_1,
sum (case when version_i > 1 then created_w_predecessors else 0 end) created_w_predecessors_v_gt1,
sum (created_w_predecessors) created_w_predecessors,
sum (case when version_i = 1 then created_wo_predecessors else 0 end) created_wo_predecessors_v_1,
sum (case when version_i > 1 then created_wo_predecessors else 0 end) created_wo_predecessors_v_gt1,

sum (created_wo_predecessors) created_wo_predecessors,
sum (active_created_prev_releases) active_created_prev_releases,
sum (active_created_this_release) active_created_this_release,
sum (active_created_next_releases) active_created_next_releases,
sum (created_this_release) created_this_release,
sum (active_updated_this_release) active_updated_this_release,
sum (active_untouched_this_release) active_untouched_this_release,
sum (active_total) active_total
from (
select
version_i,
sum (
case

  when deleted = 'Y' and last < l_prev_release
  then 1
  else 0
end) retired_prev_releases,
sum (
case
when deleted = 'Y' and last = l_prev_release then
1
else
0
end) retired_this_release,
sum (
case

  when deleted = 'Y' and last > l_prev_release
  then 1
  else 0
end) retired_next_releases,
sum (
case
when deleted = 'Y' then
1
else
0
end) retired_total,
sum (
case

  when created = l_this_release
  and exists (
    SELECT
      1
    FROM
      xref p
    WHERE
      p.ac        = x.ac
    and p.dbid    = x.dbid
    and p.created < l_this_release)
  then 1
  else 0
end) created_w_predecessors,

sum (
case
  when created = l_this_release
  and not exists (
    SELECT
      1
    FROM
      xref p
    WHERE
      p.ac        = x.ac
    and p.dbid    = x.dbid
    and p.created < l_this_release)
  then 1

  else 0
end) created_wo_predecessors,
sum (
case
  when deleted = 'N' and created < l_this_release
  then 1
  else 0
end) active_created_prev_releases,
sum (
case
  when deleted = 'N' and created = l_this_release
  then 1
  else 0

end) active_created_this_release,
sum (
case
  when deleted = 'N' and created > l_this_release
  then 1
  else 0
end) active_created_next_releases,
sum (
case
  when created = l_this_release
  then 1
  else 0
end) created_this_release,

sum (
case
  when deleted = 'N' and created != l_this_release and last = l_this_release
  then 1
  else 0
end) active_updated_this_release,
sum (
case
  when deleted = 'N' and created != l_this_release and last != l_this_release
  then 1
  else 0

end) active_untouched_this_release,

sum (
case
when deleted = 'N' then
1
else
0
end) active_total
 from xref x
where x.dbid = l_dbid
group by version_i))) q
on (s.this_release = q.this_release)
when matched then update
set

s.dbid = q.dbid,
s.prev_release = q.prev_release,
s.retired_prev_releases = q.retired_prev_releases,
s.retired_this_release = q.retired_this_release,
s.retired_next_releases = q.retired_next_releases,
s.retired_total = q.retired_total,
s.created_w_predecessors_v_1 = q.created_w_predecessors_v_1,
s.created_w_predecessors_v_gt1 = q.created_w_predecessors_v_gt1,
s.created_w_predecessors = q.created_w_predecessors,
s.created_wo_predecessors_v_1 = q.created_wo_predecessors_v_1,
s.created_wo_predecessors_v_gt1 = q.created_wo_predecessors_v_gt1,
s.created_wo_predecessors = q.created_wo_predecessors,
s.active_created_prev_releases = q.active_created_prev_releases,

s.active_created_this_release = q.active_created_this_release,
s.active_created_next_releases = q.active_created_next_releases,
s.created_this_release = q.created_this_release,
s.active_updated_this_release = q.active_updated_this_release,
s.active_untouched_this_release = q.active_untouched_this_release,
s.active_total = q.active_total
when not matched then insert
(
dbid,
this_release,
prev_release,
retired_prev_releases,
retired_this_release,

retired_next_releases,
retired_total,
created_w_predecessors_v_1,
created_w_predecessors_v_gt1,
created_w_predecessors,
created_wo_predecessors_v_1,
created_wo_predecessors_v_gt1,
created_wo_predecessors,
active_created_prev_releases,
active_created_this_release,
active_created_next_releases,
created_this_release,
active_updated_this_release,

active_untouched_this_release,
active_total)
values
(
q.dbid,
q.this_release,
q.prev_release,
q.retired_prev_releases,
q.retired_this_release,
q.retired_next_releases,
q.retired_total,
q.created_w_predecessors_v_1,
q.created_w_predecessors_v_gt1,

q.created_w_predecessors,
q.created_wo_predecessors_v_1,
q.created_wo_predecessors_v_gt1,
q.created_wo_predecessors,
q.active_created_prev_releases,
q.active_created_this_release,
q.active_created_next_releases,
q.created_this_release,
q.active_updated_this_release,
q.active_untouched_this_release,
q.active_total);
end;
/
set define on
