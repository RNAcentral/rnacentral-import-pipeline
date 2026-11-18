-- Destination tables for disease annotations (see docs/disease-ontology-plan.md).
--
--   psql -v ON_ERROR_STOP=1 -f files/schema/disease-annotations.sql "$PGDATABASE"
--
-- New, empty tables only: nothing existing is altered, so the only locks taken
-- are brief ones on ontology_terms and rnc_references while the foreign keys
-- are registered. Fails if the tables already exist (deliberately not
-- IF NOT EXISTS, so a half-applied or drifted copy is noticed).
--
-- No foreign key from urs_taxid to rnc_rna_precomputed: go_term_annotations
-- has none either, and post-release ordering does not guarantee the precomputed
-- row exists when annotations are inserted.
--
-- Undo:
--   DROP TABLE rnacen.rnc_disease_annotation_publications;
--   DROP TABLE rnacen.rnc_disease_annotations;

BEGIN;

CREATE TABLE rnacen.rnc_disease_annotations (
    disease_annotation_id bigserial PRIMARY KEY,
    urs_taxid character varying(44) NOT NULL,
    -- The source entry that made the claim.
    accession character varying(200) NOT NULL,
    -- Source database, matching rnc_database.descr.
    assigned_by character varying(50) NOT NULL,
    -- Exactly as the source wrote it.
    disease_text text NOT NULL,
    -- NULL when the text could not be mapped to a term.
    ontology_term_id character varying(15)
        REFERENCES rnacen.ontology_terms (ontology_term_id),
    match_type text NOT NULL
        CHECK (match_type IN ('override', 'label_exact', 'synonym_exact', 'unmapped')),
    CONSTRAINT un_rnc_disease_annotations
        UNIQUE (urs_taxid, accession, disease_text),
    CONSTRAINT ck_rnc_disease_annotations_mapped
        CHECK ((ontology_term_id IS NULL) = (match_type = 'unmapped'))
);

CREATE INDEX ix_rnc_disease_annotations_term
    ON rnacen.rnc_disease_annotations (ontology_term_id);
CREATE INDEX ix_rnc_disease_annotations_assigned_by
    ON rnacen.rnc_disease_annotations (assigned_by);

CREATE TABLE rnacen.rnc_disease_annotation_publications (
    disease_annotation_id bigint NOT NULL
        REFERENCES rnacen.rnc_disease_annotations (disease_annotation_id)
        ON DELETE CASCADE,
    reference_id bigint NOT NULL
        REFERENCES rnacen.rnc_references (id),
    PRIMARY KEY (disease_annotation_id, reference_id)
);

CREATE INDEX ix_rnc_disease_annotation_publications_ref
    ON rnacen.rnc_disease_annotation_publications (reference_id);

COMMIT;
