DROP TABLE IF EXISTS load_rnc_accessions;
CREATE TABLE load_rnc_accessions (
    accession character varying(200) NOT NULL,
    parent_ac character varying(200) NULL,
    seq_version bigint NULL,
    feature_start bigint NULL,
    feature_end bigint NULL,
    feature_name character varying(80) NULL,
    ordinal bigint NULL,
    division character varying(6) NULL,
    keywords character varying(200) NULL,
    description character varying(500) NULL,
    species character varying(300) NULL,
    organelle character varying(200) NULL,
    classification character varying(1000) NULL,
    project character varying(100) NULL,
    is_composite character varying(2) NULL,
    non_coding_id character varying(200) NULL,
    database character varying(40) NULL,
    external_id character varying(300) NULL,
    optional_id character varying(200) NULL,
    common_name character varying(200) NULL,
    allele character varying(100) NULL,
    anticodon character varying(200) NULL,
    chromosome character varying(200) NULL,
    experiment text NULL,
    function character varying(4000) NULL,
    gene character varying(200) NULL,
    gene_synonym character varying(800) NULL,
    inference character varying(600) NULL,
    locus_tag character varying(100) NULL,
    map character varying(400) NULL,
    mol_type character varying(100) NULL,
    ncrna_class character varying(100) NULL,
    note text NULL,
    old_locus_tag character varying(200) NULL,
    operon character varying(100) NULL,
    product character varying(600) NULL,
    pseudogene character varying(100) NULL,
    standard_name character varying(200) NULL,
    db_xref text NULL,
    so_term character varying(15) NULL,
    url text NULL
);

DROP TABLE IF EXISTS load_rnc_references;
CREATE TABLE load_rnc_references (
    md5 character varying(64) NOT NULL,
    accession character varying(200) NULL,
    authors text NULL,
    location character varying(4000) NULL,
    title character varying(4000) NULL,
    pmid character varying(20) NULL,
    doi text NULL
  );

DROP TABLE IF EXISTS load_rnacentral;
CREATE TABLE
  load_rnacentral (
    crc64 character varying(16) NULL,
    len integer NULL,
    seq_short character varying(4000) NULL,
    seq_long text NULL,
    database character varying(40) NULL,
    ac character varying(200) NULL,
    optional_id character varying(100) NULL,
    version integer NULL,
    taxid bigint NULL,
    md5 character varying(32) NULL
  );

DROP TABLE IF EXISTS load_rnacentral_all;
CREATE TABLE
  load_rnacentral_all (
    crc64 character varying(16) NULL,
    len integer NULL,
    seq_short character varying(4000) NULL,
    seq_long text NULL,
    database character varying(40) NULL,
    ac character varying(200) NULL,
    optional_id character varying(100) NULL,
    version integer NULL,
    taxid bigint NULL,
    md5 character varying(32) NULL
  );

DROP TABLE IF EXISTS load_assemblies;
CREATE TABLE load_assemblies (
  assembly_id varchar(255) NOT NULL,
  assembly_full_name varchar(255) NOT NULL,
  gca_accession varchar(20) NULL,
  assembly_ucsc varchar(100) NULL,
  common_name varchar(255),
  taxid int4 NOT NULL,
  ensembl_url varchar(100) NULL,
  division varchar(20) NULL,
  blat_mapping int4 NULL,
  example_chromosome varchar(20) NULL,
  example_end int4 NULL,
  example_start int4 NULL,
  subdomain varchar(100) NOT NULL
);

DROP TABLE IF EXISTS load_compara;
CREATE TABLE load_compara (
  homology_group text not null,
  ensembl_transcript text not null
);

DROP TABLE IF EXISTS load_coordinate_info;
CREATE TABLE load_coordinate_info (
  chromosome text NOT NULL,
  coordinate_system text NOT NULL,
  assembly_id text,
  is_reference bool,
  karyotype_rank int
);

DROP TABLE IF EXISTS load_ensembl_analysis_status;
CREATE TABLE load_ensembl_analysis_status (
    database_name text,
    task_name text
);

DROP TABLE IF EXISTS load_genome_mapping_attempted;
CREATE TABLE load_genome_mapping_attempted (
  urs_taxid text not null,
  assembly_id text not null,
  last_run timestamp not null default CURRENT_TIMESTAMP
);

DROP TABLE IF EXISTS load_assemblies;
CREATE TABLE load_assemblies (
  assembly_id varchar(255) NOT NULL,
  assembly_full_name varchar(255) NOT NULL,
  gca_accession varchar(20) NULL,
  assembly_ucsc varchar(100) NULL,
  common_name varchar(255),
  taxid int4 NOT NULL,
  ensembl_url varchar(100) NULL,
  division varchar(20) NULL,
  blat_mapping int4 NULL,
  example_chromosome varchar(20) NULL,
  example_end int4 NULL,
  example_start int4 NULL,
  subdomain varchar(100) NOT NULL
);

DROP TABLE IF EXISTS load_rfam_model_hits;
CREATE TABLE load_rfam_model_hits (
  sequence_start integer NOT NULL,
  sequence_stop integer NOT NULL,
  sequence_completeness double precision,
  model_start integer NOT NULL,
  model_stop integer NOT NULL,
  model_completeness double precision,
  overlap character varying(30) COLLATE pg_catalog."default" NOT NULL,
  e_value double precision NOT NULL,
  score double precision NOT NULL,
  rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
  upi character varying(13) COLLATE pg_catalog."default" NOT NULL
);

DROP TABLE IF EXISTS load_rfam_model_hits;
CREATE TABLE load_rfam_model_hits (
  urs text PRIMARY KEY REFERENCES rna(upi),
  qa_analysis text NOT NULL,
  last_run timestamp NOT NULL
);

DROP TABLE IF EXISTS load_rnc_sequence_features;
CREATE TABLE load_rnc_sequence_features (
    accession varchar(100) NOT NULL,
    taxid int not null,
    start int not null,
    stop int not null,
    feature_name varchar(50),
    metadata jsonb
);

DROP TABLE IF EXISTS load_go_term_annotations;
CREATE TABLE load_go_term_annotations (
    rna_id varchar(50),
    qualifier text,
    ontology_term_id varchar(15),
    evidence_code varchar(15),
    assigned_by varchar(50),
    extensions jsonb
);

DROP TABLE IF EXISTS load_go_term_publication_map;
CREATE TABLE load_go_term_publication_map (
    rna_id varchar(50),
    qualifier text,
    assigned_by varchar(50),
    ontology_term_id varchar(15),
    evidence_code varchar(15),
    pubmed_id text
);

DROP TABLE IF EXISTS load_ref_pubmed;
CREATE TABLE load_ref_pubmed (
    ref_pubmed_id int,
    authors text,
    location text,
    title text,
    doi text
);

DROP TABLE IF EXISTS load_interactions;
CREATE TABLE load_interactions (
  intact_id text NOT NULL,
  urs_taxid text NOT NULL,
  interacting_id text NOT NULL,
  names JSONB NOT NULL,
  taxid int NOT NULL
);

DROP TABLE IF EXISTS load_karyotypes;
CREATE TABLE load_karyotypes (
	assembly_id varchar(255) NOT NULL,
    karyotype text
);

DROP TABLE IF EXISTS load_rnc_coordinates;
CREATE TABLE load_rnc_coordinates (
    accession varchar(200) NULL,
    local_start int8 NULL,
    local_end int8 NULL,
    chromosome varchar(100) NULL,
    strand int8 NULL,
    assembly_id varchar(255) NULL
);

DROP TABLE IF EXISTS load_ontology_terms;
CREATE TABLE load_ontology_terms (
  ontology_term_id varchar(15),
  ontology varchar(5),
  name text,
  definition text
);

DROP TABLE IF EXISTS load_protein_info;
CREATE TABLE load_protein_info (
  protein_accession text NOT NULL,
  description text,
  label text,
  synonyms text[]
);

DROP TABLE IF EXISTS load_rnc_sequence_regions;
CREATE TABLE load_rnc_sequence_regions (
    accession text,
    urs_taxid text,
    region_name text not null,
    chromosome text,
    strand int4,
    exon_start int4,
    exon_stop int4,
    assembly_id varchar(255),
    exon_count int
);

DROP TABLE IF EXISTS load_rnc_related_sequences;
CREATE TABLE load_rnc_related_sequences (
  source_accession varchar(100) NOT NULL,
  source_urs_taxid text,
  target_accession varchar(100) NOT NULL,
  relationship_type text NOT NULL,
  methods text[]
);

DROP TABLE IF EXISTS load_rfam_clans;
CREATE TABLE load_rfam_clans (
    rfam_clan_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    name character varying(40) COLLATE pg_catalog."default" NOT NULL,
    description character varying(1000) COLLATE pg_catalog."default" NOT NULL,
    family_count integer NOT NULL
);

DROP TABLE IF EXISTS load_rfam_models;
CREATE TABLE load_rfam_models (
    rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
    long_name character varying(200) COLLATE pg_catalog."default" NOT NULL,
    description character varying(2000) COLLATE pg_catalog."default",
    seed_count integer NOT NULL,
    full_count integer NOT NULL,
    length integer NOT NULL,
    is_suppressed boolean NOT NULL,
    rfam_clan_id character varying(20) COLLATE pg_catalog."default",
    domain character varying(50) COLLATE pg_catalog."default",
    rna_type character varying(250) COLLATE pg_catalog."default",
    short_name character varying(50) COLLATE pg_catalog."default",
    rfam_rna_type character varying(250) COLLATE pg_catalog."default",
    so_rna_type text
);

DROP TABLE IF EXISTS load_rfam_go_terms;
CREATE TABLE load_rfam_go_terms (
    ontology_term_id character varying(10) COLLATE pg_catalog."default" NOT NULL,
    rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL
);

DROP TABLE IF EXISTS load_rnc_secondary_structure;
CREATE TABLE load_rnc_secondary_structure (
    rnc_accession_id varchar(100),
    secondary_structure text,
    md5 varchar(32)
);

DROP TABLE IF EXISTS load_taxonomy;
CREATE TABLE load_taxonomy (
    taxid int,
    name text,
    lineage text,
    aliases json,
    replaced_by int
);

DROP TABLE IF EXISTS load_overlaps;
CREATE TABLE load_overlaps (
    upi_taxid text,
    status text,
    result text,
    assembly_id text
);

DROP TABLE IF EXISTS load_qa_status;
CREATE TABLE load_qa_status (
  rna_id varchar(44) NOT NULL,
  upi varchar(26) NOT NULL,
  taxid int8 NOT NULL,
  has_issue bool,
  incomplete_sequence bool,
  possible_contamination bool,
  missing_rfam_match bool,
  messages jsonb
);

DROP TABLE IF EXISTS load_qa_rfam_attempted;
CREATE TABLE load_qa_rfam_attempted (
  urs text NOT NULL,
  model_source text NOT NULL,
  source_version text NOT NULL,
  last_run timestamp default CURRENT_TIMESTAMP
);

DROP TABLE IF EXISTS load_rfam_model_hits;
CREATE TABLE load_rfam_model_hits (
  sequence_start integer NOT NULL,
  sequence_stop integer NOT NULL,
  sequence_completeness double precision,
  model_start integer NOT NULL,
  model_stop integer NOT NULL,
  model_completeness double precision,
  overlap character varying(30) COLLATE pg_catalog."default" NOT NULL,
  e_value double precision NOT NULL,
  score double precision NOT NULL,
  rfam_model_id character varying(20) COLLATE pg_catalog."default" NOT NULL,
  upi character varying(13) COLLATE pg_catalog."default" NOT NULL
);

DROP TABLE IF EXISTS load_rnc_text_mining;
CREATE TABLE load_rnc_text_mining (
    pattern_group text,
    pattern text,
    matching_word text,
    sentence text,
    md5 text,
    authors text,
    location text,
    title text,
    pmid text,
    doi text
);

DROP TABLE IF EXISTS load_secondary_layout_models;
CREATE TABLE load_secondary_layout_models (
    model_name text NOT NULL,
    taxid int NOT NULL,
    rna_type text NOT NULL,
    so_term text NOT NULL,
    cell_location text NOT NULL,
    model_source text not null,
    model_length int not null
);

DROP TABLE IF EXISTS load_secondary_layout;
CREATE TABLE load_secondary_layout (
    urs text NOT NULL,
    secondary_structure text NOT NULL,
    model text NOT NULL,
    overlap_count int not null,
    basepair_count int not null,
    model_start int,
    model_stop int,
    model_coverage float,
    sequence_start int,
    sequence_stop int,
    sequence_coverage float
);

DROP TABLE IF EXISTS load_ensembl_pseudogenes;
create table load_ensembl_pseudogenes (
    gene text not null,
    region_name text not null,
    chromosome text not null,
    strand int4 not null,
    exon_start int4 not null,
    exon_stop int4 not null,
    assembly_id varchar(255) not null,
    exon_count int not null
);
