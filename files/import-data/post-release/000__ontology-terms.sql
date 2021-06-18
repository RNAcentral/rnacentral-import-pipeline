BEGIN TRANSACTION;

insert into ontology_terms (
  ontology_term_id,
  ontology,
  name,
  definition
) (
select distinct
  ontology_term_id,
  ontology,
  name,
  definition
from load_ontology_terms
)
ON CONFLICT (ontology_term_id) DO UPDATE SET
  ontology_term_id = excluded.ontology_term_id,
  ontology = excluded.ontology,
  name = excluded.name,
  definition = excluded.definition
;

drop table load_ontology_terms;

COMMIT;
