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
    products: Option<String>,
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

impl FromIterator<RawAccession> for AccessionVec {
    fn from_iter<I: IntoIterator<Item = RawAccession>>(iter: I) -> Self {
        let mut a = AccessionVec::default();

        for i in iter {
            a.species.extend(i.species.into_iter().filter(|s| !s.is_empty()));
            a.tax_strings.extend(i.tax_strings.into_iter().filter(|t| !t.is_empty()));
            a.organelles.extend(i.organelles.clone().filter(|t| !t.is_empty()));
            a.functions.extend(i.functions.into_iter().filter(|t| !t.is_empty()));
            a.genes.extend(i.genes.into_iter().filter(|t| !t.is_empty()));
            a.common_name.extend(i.common_name.into_iter().filter(|t| !t.is_empty()));
            a.notes.extend(i.notes.into_iter().filter(|t| !t.is_empty()));
            a.locus_tags.extend(i.locus_tags.into_iter().filter(|t| !t.is_empty()));
            a.standard_names.extend(i.standard_names.into_iter().filter(|t| !t.is_empty()));
            a.products.extend(i.products.into_iter().filter(|t| !t.is_empty()));

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
            if i.authors.is_some() {
                let authors = i.authors.unwrap();
                let authors = authors.split(", ").filter(|a| !a.is_empty()).map(|s| s.to_string());
                value.authors.extend(authors);
            }

            value.journals.extend(i.journal.into_iter());
            value.pub_titles.extend(i.pub_title.into_iter().filter(|t| !t.is_empty()));
            value.pub_ids.extend(i.pub_id.into_iter());
            value.pubmed_ids.extend(i.pubmed_id.into_iter().filter(|t| !t.is_empty()));
            value.dois.extend(i.doi.into_iter().filter(|t| !t.is_empty()));
        }

        value
    }
}
