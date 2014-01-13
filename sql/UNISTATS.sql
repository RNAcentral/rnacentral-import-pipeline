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

create or replace PACKAGE                                 UNISTATS
-- authid current_user
is

-- gather ONLY SUBPARTITION level stats
  PROCEDURE gather_subpartition_stats (
     p_tabname varchar2,
     p_subpartname varchar2,
     p_ownname varchar2 default user
  );

-- gather ONLY PARTITION level stats
  PROCEDURE gather_partition_stats (

     p_tabname varchar2,
     p_partname varchar2,
     p_ownname varchar2 default user
  );

-- gather ONLY GLOBAL level stats
  PROCEDURE gather_global_table_stats (
     p_tabname varchar2,
     p_ownname varchar2 default user
  );

-- gather ALL levels stats
  PROCEDURE gather_table_stats (

     p_tabname varchar2,
     p_ownname varchar2 default user
  );

END UNISTATS;
/
create or replace PACKAGE BODY                                         UNISTATS AS


  procedure gather_table_stats_internal (
     p_ownname varchar2,
     p_tabname varchar2,
     p_partname varchar2 default null,

     p_estimate_percent number default dbms_stats.AUTO_SAMPLE_SIZE,
     p_block_sample boolean default FALSE,
     p_method_opt varchar2 default dbms_stats.DEFAULT_METHOD_OPT,
     p_degree number default dbms_stats.AUTO_DEGREE,
     p_granularity varchar2 default  dbms_stats.DEFAULT_GRANULARITY,
     p_cascade boolean default dbms_stats.AUTO_CASCADE,
     p_stattab varchar2 default null,
     p_statid varchar2 default null,
     p_statown varchar2 default null,
     p_no_invalidate boolean default dbms_stats.DEFAULT_NO_INVALIDATE,
     p_stattype varchar2 default 'DATA',
     p_force boolean default FALSE
  ) AS

  BEGIN
    dbms_stats.gather_table_stats (
      ownname => p_ownname,
      tabname => p_tabname,
      partname => p_partname,
      estimate_percent => p_estimate_percent,
      block_sample => p_block_sample,
      method_opt => p_method_opt,
      degree => p_degree,
      granularity => p_granularity,
      cascade => p_cascade,
      stattab => p_stattab,
      statid => p_statid,

      statown => p_statown,
      no_invalidate => p_no_invalidate,
      stattype => p_stattype,
      force => p_force
    );

  END gather_table_stats_internal;

  PROCEDURE gather_subpartition_stats (
     p_tabname varchar2,
     p_subpartname varchar2,
     p_ownname varchar2 default user
  )

  AS
  begin
    gather_table_stats_internal (
      p_ownname => p_ownname,
      p_tabname => p_tabname,
      p_partname => p_subpartname,
      p_estimate_percent => dbms_stats.AUTO_SAMPLE_SIZE,
      p_method_opt => 'FOR ALL COLUMNS SIZE AUTO',
      p_degree => DBMS_STATS.AUTO_DEGREE,
      p_granularity => 'SUBPARTITION',
      p_cascade => TRUE
    );
  end gather_subpartition_stats;


  PROCEDURE gather_partition_stats (
     p_tabname varchar2,
     p_partname varchar2,
     p_ownname varchar2 default user
  )
  AS
  begin
    gather_table_stats_internal (
      p_ownname => p_ownname,
      p_tabname => p_tabname,
      p_partname => p_partname,
      p_estimate_percent => dbms_stats.AUTO_SAMPLE_SIZE,

      p_method_opt => 'FOR ALL COLUMNS SIZE AUTO',
      p_degree => DBMS_STATS.AUTO_DEGREE,
      p_granularity => 'APPROX_GLOBAL AND PARTITION',
      p_cascade => TRUE
    );
  end gather_partition_stats;

  PROCEDURE gather_global_table_stats (
     p_tabname varchar2,
     p_ownname varchar2 default user
  )
  AS
  begin

    gather_table_stats_internal (
      p_ownname => p_ownname,
      p_tabname => p_tabname,
      p_estimate_percent => dbms_stats.AUTO_SAMPLE_SIZE,
      p_method_opt => 'FOR ALL COLUMNS SIZE AUTO',
      p_degree => DBMS_STATS.AUTO_DEGREE,
      p_granularity => 'GLOBAL',
      p_cascade => TRUE

    );
  end gather_global_table_stats;

  PROCEDURE gather_table_stats (

     p_tabname varchar2,
     p_ownname varchar2 default user
  )
  AS
  begin
    gather_table_stats_internal (
      p_ownname => p_ownname,
      p_tabname => p_tabname,
      p_estimate_percent => dbms_stats.AUTO_SAMPLE_SIZE,
      p_method_opt => 'FOR ALL COLUMNS SIZE AUTO',
      p_degree => DBMS_STATS.AUTO_DEGREE,
      p_granularity => 'ALL',
      p_cascade => TRUE

    );
  end gather_table_stats;

END UNISTATS;
/
set define on
