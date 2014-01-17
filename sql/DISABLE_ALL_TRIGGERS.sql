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

create or replace PROCEDURE        "DISABLE_ALL_TRIGGERS"
IS
/***********************************************************************
/* Package Name:  DISABLE_ALL_TRIGGERS
/* Author:        Steven Rosanoff
/* Date Created:  26/10/2012
/* Description:   This package enables all triggers belonging to UNIPARC.
/* Automated:     NO
/**********************************************************************/

    CURSOR trigger_list
    IS
      SELECT trigger_name

        FROM sys.all_triggers
       WHERE owner = 'UNIPARC';

  BEGIN

      FOR trig IN trigger_list
      LOOP
        EXECUTE IMMEDIATE
        'ALTER TRIGGER ' || trig.trigger_name || ' DISABLE';
      END LOOP;

  END DISABLE_ALL_TRIGGERS;
/
set define on
