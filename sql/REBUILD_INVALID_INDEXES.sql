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

create or replace PROCEDURE        "REBUILD_INVALID_INDEXES"
IS
/***********************************************************************
/* Package Name:  REBUILD_INVALID_INDEXES
/* Author:        Steven Rosanoff
/* Date Created:  13/12/2012
/* Description:   This package rebuilds any indexes marked as unusable.
/* Automated:     NO. This procedure is called when the new release is
/*                prepared. The script that calls this procedure is:
/**********************************************************************/

    CURSOR index_list
    IS

      SELECT INDEX_NAME
      FROM sys.ALL_INDEXES
      WHERE status= 'UNUSABLE';

  BEGIN

      FOR ind IN index_list
      LOOP
        EXECUTE IMMEDIATE
        'ALTER INDEX ' || ind.index_name || ' REBUILD TABLESPACE UNIIND';
      END LOOP;

  END REBUILD_INVALID_INDEXES;
/
set define on
