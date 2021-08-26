use std::{
    collections::HashSet,
    iter::FromIterator,
};

use serde::{
    Deserialize,
    Serialize,
};

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, Clone)]
pub struct RawAccession {
    pub id: usize,
    pub urs_taxid: String,
    accession: String,
    common_name: Option<String>,
    database: String,
    external_id: String,
    function: Option<String>,
    gene_synonyms: Option<String>,
    gene: Option<String>,
    locus_tag: Option<String>,
    non_coding_id: Option<String>,
    notes: Option<String>,
    optional_id: Option<String>,
    organelle: Option<String>,
    parent_accession: Option<String>,
    product: Option<String>,
    species: Option<String>,
    standard_name: Option<String>,
    tax_string: Option<String>,
    authors: Option<String>,
    journal: Option<String>,
    pub_title: Option<String>,
    pub_id: Option<u64>,
    pubmed_id: Option<String>,
    doi: Option<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, Clone)]
pub struct CrossReference {
    name: String,
    external_id: String,
    optional_id: Option<String>,
    accession: String,
    non_coding_id: Option<String>,
    parent_accession: Option<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct AccessionVec {
    species: HashSet<String>,
    organelles: HashSet<String>,
    tax_strings: HashSet<String>,
    functions: HashSet<String>,
    genes: HashSet<String>,
    gene_synonyms: HashSet<String>,
    common_name: HashSet<String>,
    notes: HashSet<String>,
    locus_tags: HashSet<String>,
    standard_names: HashSet<String>,
    products: HashSet<String>,
}

#[derive(Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReferenceVec {
    authors: HashSet<String>,
    journals: HashSet<String>,
    pub_titles: HashSet<String>,
    pub_ids: HashSet<u64>,
    pubmed_ids: HashSet<String>,
    dois: HashSet<String>,
}

impl Default for ReferenceVec {
    fn default() -> Self {
        Self {
            authors: HashSet::new(),
            journals: HashSet::new(),
            pub_titles: HashSet::new(),
            pub_ids: HashSet::new(),
            pubmed_ids: HashSet::new(),
            dois: HashSet::new(),
        }
    }
}

impl Default for AccessionVec {
    fn default() -> Self {
        Self {
            species: HashSet::new(),
            organelles: HashSet::new(),
            tax_strings: HashSet::new(),
            functions: HashSet::new(),
            genes: HashSet::new(),
            gene_synonyms: HashSet::new(),
            common_name: HashSet::new(),
            notes: HashSet::new(),
            locus_tags: HashSet::new(),
            standard_names: HashSet::new(),
            products: HashSet::new(),
        }
    }
}

impl From<RawAccession> for CrossReference {
    fn from(raw: RawAccession) -> Self {
        Self {
            name: raw.database,
            external_id: raw.external_id,
            optional_id: raw.optional_id,
            accession: raw.accession,
            non_coding_id: raw.non_coding_id,
            parent_accession: raw.parent_accession,
        }
    }
}

fn insert(current: &mut HashSet<String>, raw: Option<String>) {
    match raw {
        None => (),
        Some(s) => match s.is_empty() {
            false => {
                current.insert(s);
            },
            true => (),
        },
    }
}

impl FromIterator<RawAccession> for AccessionVec {
    fn from_iter<I: IntoIterator<Item = RawAccession>>(iter: I) -> Self {
        let mut a = AccessionVec::default();

        for i in iter {
            insert(&mut a.species, i.species);
            insert(&mut a.tax_strings, i.tax_string);
            insert(&mut a.functions, i.function);
            insert(&mut a.organelles, i.organelle.map(|s| s.to_lowercase()));
            insert(&mut a.genes, i.gene);
            insert(&mut a.common_name, i.common_name.map(|s| s.to_lowercase()));
            insert(&mut a.notes, i.notes);
            insert(&mut a.locus_tags, i.locus_tag);
            insert(&mut a.standard_names, i.standard_name);
            insert(&mut a.products, i.product);

            for synonyms in i.gene_synonyms.into_iter() {
                if synonyms.contains(";") {
                    a.gene_synonyms
                        .extend(synonyms.split(";").map(|p| p.trim()).map(str::to_string))
                } else if synonyms.contains(",") {
                    a.gene_synonyms
                        .extend(synonyms.split(",").map(|p| p.trim()).map(str::to_string))
                } else {
                    a.gene_synonyms.insert(synonyms);
                }
            }
        }

        a
    }
}

impl FromIterator<RawAccession> for ReferenceVec {
    fn from_iter<I: IntoIterator<Item = RawAccession>>(iter: I) -> Self {
        let mut value = ReferenceVec::default();

        for i in iter {
            match i.authors {
                None => (),
                Some(authors) => {
                    let authors =
                        authors.split(", ").filter(|a| !a.is_empty()).map(|s| s.to_string());
                    value.authors.extend(authors);
                },
            }

            insert(&mut value.journals, i.journal);
            insert(&mut value.pub_titles, i.pub_title);
            insert(&mut value.pubmed_ids, i.pubmed_id);
            insert(&mut value.dois, i.doi);
            value.pub_ids.extend(i.pub_id.into_iter());

            value.pub_titles.remove("None");
        }

        value
    }
}

impl CrossReference {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn external_id(&self) -> &str {
        &self.external_id
    }
}
