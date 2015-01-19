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

create or replace PACKAGE RNC_HEALTHCHECKS AS

  procedure run_healthchecks;

END RNC_HEALTHCHECKS;
/
create or replace PACKAGE BODY RNC_HEALTHCHECKS AS


  /*
  * All UPIs must be unique.
  */
  PROCEDURE check_unique_upis AS

    v_distinct_upi NUMBER;
    v_distinct_md5 NUMBER;
    v_all_rna_seq  NUMBER;
  BEGIN

    SELECT count(DISTINCT upi) INTO v_distinct_upi FROM rna;
    SELECT count(DISTINCT md5) INTO v_distinct_md5 FROM rna;
    SELECT count(*)            INTO v_all_rna_seq  FROM rna;

    IF v_distinct_upi != v_distinct_md5 AND
       v_distinct_md5 != v_all_rna_seq THEN
      DBMS_OUTPUT.put_line('not ok ... check_unique_upis');
    ELSE

      DBMS_OUTPUT.put_line('ok ... check_unique_upis, found '
                            || v_all_rna_seq || ' RNA sequences');
    END IF;

  END check_unique_upis;


  /*
  * Taxids should not be negative. This problem was encountered once.
  */
  PROCEDURE check_negative_taxid AS
    v_xref NUMBER;
    v_staging NUMBER;

  BEGIN
    SELECT count(*) INTO v_xref FROM xref WHERE taxid < 0;
    SELECT count(*) into v_staging FROM load_rnacentral_all WHERE taxid < 0;

    IF v_xref != v_staging OR v_xref != 0 THEN
      DBMS_OUTPUT.put_line('not ok ... check_negative_taxid');
    ELSE
      DBMS_OUTPUT.put_line('ok ... check_negative_taxid');
    end if;

  end check_negative_taxid;



  /*
  * Sometimes the staging table contains duplicate accession numbers.
  * The update should not start if this happens.
  */
  procedure check_unique_ac_staging
  AS
    v_ac_all NUMBER;
    v_ac_distinct NUMBER;
  BEGIN
    SELECT count(*) into v_ac_all FROM load_rnacentral_all;
    SELECT count(DISTINCT ac) into v_ac_distinct FROM load_rnacentral_all;

    IF v_ac_all != v_ac_distinct THEN

      DBMS_OUTPUT.put_line('not ok ... check_unique_ac_staging');
    ELSE
      DBMS_OUTPUT.put_line('ok ... check_unique_ac_staging');
    end if;

  end check_unique_ac_staging;

  /*
  * Make sure that all accessions in xref are unique.
  */
  procedure check_unique_ac_xref
  AS
    v_ac_all NUMBER;

    v_ac_distinct NUMBER;
  BEGIN
    SELECT count(*) INTO v_ac_all FROM xref WHERE deleted='N';
    SELECT count(DISTINCT ac) into v_ac_distinct FROM xref where deleted='N';

    IF v_ac_all != v_ac_distinct THEN
      DBMS_OUTPUT.put_line('not ok ... check_unique_ac_xref');
    ELSE
      DBMS_OUTPUT.put_line('ok ... check_unique_ac_xref');
    end if;

  END check_unique_ac_xref;



  PROCEDURE check_rnas_without_xrefs
  AS
    v_count NUMBER;
  BEGIN
    SELECT count(*) INTO v_count FROM xref x, rna r WHERE x.upi (+) = r.upi AND x.upi IS NULL;

    IF v_count != 0 THEN
      DBMS_OUTPUT.put_line('not ok ... check_rnas_without_xrefs');
    else
      DBMS_OUTPUT.put_line('ok ... check_rnas_without_xrefs');
    END IF;


  END check_rnas_without_xrefs;


  PROCEDURE check_xrefs_without_lit_refs
  AS
    v_count NUMBER;
  BEGIN
    select count(distinct upi) INTO v_count
    from xref t1, rnc_reference_map t2
    where
    t2.ACCESSION (+) = t1.ac and
    t2.accession is null and
    deleted = 'N';

    -- some entries don't have literature refs in ENA
    IF v_count > 100 THEN
      DBMS_OUTPUT.put_line('not ok ... check_xrefs_without_literature_refs');
    else
      DBMS_OUTPUT.put_line('ok ... check_xrefs_without_literature_refs');
    END IF;

  END check_xrefs_without_lit_refs;

  PROCEDURE check_xrefs_without_ac_data
  AS
    v_count NUMBER;
  BEGIN
    SELECT count(*) INTO v_count
    FROM xref t1, rnc_accessions t2
    WHERE
      t2.ACCESSION (+) = t1.ac AND
      t2.accession is null;

    IF v_count > 0 THEN
      DBMS_OUTPUT.put_line('not ok ... check_xrefs_without_ac_data');
    else
      DBMS_OUTPUT.put_line('ok ... check_xrefs_without_ac_data');
    END IF;

  END check_xrefs_without_ac_data;

  /*
  * Auxiliary procedure to report healthcheck results.
  */
  PROCEDURE report_results (
    p_in_count IN NUMBER,
    p_in_msg IN STRING
  )
  AS
  BEGIN
    IF p_in_count > 0 THEN
      DBMS_OUTPUT.put_line('not ok ... '||p_in_msg);
    ELSE
      DBMS_OUTPUT.put_line('ok ... '||p_in_msg);
    END IF;

  END report_results;


  /*
  * Verify external and optional ids in rnc_accessions.
  * Add more checks as more expert databases are imported into RNAcentral.
  */
  PROCEDURE check_expert_db_accessions
  AS
    v_count NUMBER;
  BEGIN

    -- RFAM
    SELECT count(*) INTO v_count FROM rnc_accessions WHERE database = 'RFAM' AND (external_id IS NULL OR optional_id IS NULL);
    report_results(v_count, 'RFAM expert accessions');

    -- SRPDB, only external_id is provided
    SELECT count(*) INTO v_count FROM rnc_accessions WHERE database = 'SRPDB' AND external_id IS null;
    report_results(v_count, 'SRPDB expert accessions');

    -- MIRBASE
    SELECT count(*) INTO v_count FROM rnc_accessions WHERE database = 'MIRBASE' AND (external_id IS NULL OR optional_id IS NULL);
    report_results(v_count, 'MIRBASE expert accessions');

    -- VEGA
    SELECT count(*) INTO v_count FROM rnc_accessions WHERE database = 'VEGA' AND (external_id IS NULL OR optional_id IS NULL);
    report_results(v_count, 'VEGA expert accessions');

    -- tmRNA Website, not all entries have both external and optional ids
    SELECT count(*) into v_count FROM rnc_accessions WHERE database = 'TMRNA_WEB' AND external_id IS null;
    report_results(v_count, 'TMRNA_WEB expert accessions');

    -- gtRNAdb, only external ids
    SELECT count(*) INTO v_count FROM rnc_accessions WHERE database = 'GTRNADB' AND (external_id IS NULL);
    report_results(v_count, 'GTRNADB expert accessions');

    -- LNCRNADB
    SELECT count(*) INTO v_count FROM rnc_accessions WHERE database = 'LNCRNADB' AND (external_id IS NULL OR optional_id IS NULL);
    report_results(v_count, 'LNCRNADB expert accessions');

  END check_expert_db_accessions;


  /*
  * After full update the id column is set to NULL
  * because the xrefs are imported using Partition Exchange Loading so no triggers are fired,
  * but this prevents django models to work correctly, so the ids must be populated before
  * the data are used by the web code.
  */
  PROCEDURE check_xref_id_not_null
  AS
    v_count NUMBER;
  BEGIN
    SELECT count(*) INTO v_count FROM xref WHERE id IS NULL;
    report_results(v_count, 'check_xref_id_not_null');
  END check_xref_id_not_null;

  /*
  * Main entry point.
  */
  PROCEDURE run_healthchecks AS
  BEGIN

    check_unique_upis;
    check_negative_taxid;
    check_unique_ac_staging;
    check_unique_ac_xref;
    check_rnas_without_xrefs;
    check_xrefs_without_lit_refs;
    check_xrefs_without_ac_data;
    check_expert_db_accessions;
    check_xref_id_not_null;

  END run_healthchecks;



END RNC_HEALTHCHECKS;
/
set define on
