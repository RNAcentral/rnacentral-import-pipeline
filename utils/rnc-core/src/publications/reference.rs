use thiserror::Error;

use md5::{
    Digest,
    Md5,
};

use crate::publications::external_reference::ExternalReference;

#[derive(Error, Debug)]
pub enum ReferenceBuildError {
    #[error("Reference must have a title")]
    NoTitle,
}

#[derive(Error, Debug)]
pub enum AuthorBuildingError {
    #[error("Author must have a name")]
    NoName,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Author(String, String);

#[derive(Debug, PartialEq, Eq)]
pub struct RefExternal {
    pmid: Option<String>,
    pmcid: Option<String>,
    doi: Option<String>,
}

#[derive(Debug, PartialEq, Eq)]
pub struct Reference {
    title: String,
    authors: Vec<Author>,
    journal: String,
    year: String,
    external: RefExternal,
}

pub struct AuthorBuilder {
    first: Option<String>,
    last: Option<String>,
}

pub struct ReferenceBuilder {
    title: Option<String>,
    authors: Vec<Author>,
    journal: Option<String>,
    year: Option<String>,
    pmid: Option<String>,
    doi: Option<String>,
    pmcid: Option<String>,
}

impl Reference {
    pub fn builder() -> ReferenceBuilder {
        ReferenceBuilder::new()
    }

    pub fn external_ids(&self) -> Vec<ExternalReference> {
        let mut ids = Vec::with_capacity(3);
        let raw = vec![self.pmid(), self.doi(), self.pmcid()];
        for external in raw {
            if external.is_some() {
                ids.push(external.unwrap())
            }
        }
        ids
    }

    pub fn location(&self) -> String {
        "".to_string()
    }

    pub fn authors(&self) -> String {
        "".to_string()
    }

    pub fn title(&self) -> &String {
        &self.title
    }

    pub fn year(&self) -> &String {
        &self.year
    }

    pub fn pmid(&self) -> Option<ExternalReference> {
        self.external.pmid.as_ref().map(|s| ExternalReference::pmid(s.to_string()))
    }

    pub fn pmcid(&self) -> Option<ExternalReference> {
        self.external.pmcid.as_ref().map(|s| ExternalReference::pmcid(s.to_string()))
    }

    pub fn doi(&self) -> Option<ExternalReference> {
        self.external.doi.as_ref().map(|s| ExternalReference::doi(s.to_string()))
    }

    pub fn md5(&self) -> String {
        let mut hasher = Md5::new();
        hasher.update(self.authors().as_bytes());
        hasher.update(self.location().as_bytes());
        hasher.update(self.title().as_bytes());
        format!("{:x}", hasher.finalize())
    }
}

impl AuthorBuilder {
    pub fn new() -> Self {
        Self {
            last: None,
            first: None,
        }
    }

    pub fn set_last_name(&mut self, last: String) {
        self.last = Some(last);
    }

    pub fn set_first_name(&mut self, first: String) {
        self.first = Some(first);
    }

    pub fn build(self) -> Result<Author, AuthorBuildingError> {
        if self.first.is_none() && self.last.is_none() {
            return Err(AuthorBuildingError::NoName);
        }

        Ok(Author(self.first.unwrap_or("".to_string()), self.last.unwrap_or("".to_string())))
    }
}

impl From<(String, String)> for Author {
    fn from(given: (String, String)) -> Self {
        let (first, last) = given;
        Author(first, last)
    }
}

impl From<(&str, &str)> for Author {
    fn from(given: (&str, &str)) -> Self {
        let (first, last) = given;
        Author(first.to_string(), last.to_string())
    }
}

impl ReferenceBuilder {
    pub fn new() -> Self {
        Self {
            title: None,
            authors: Vec::new(),
            journal: None,
            year: None,
            pmid: None,
            pmcid: None,
            doi: None,
        }
    }

    pub fn set_title(&mut self, title: String) {
        self.title = Some(title);
    }

    pub fn add_author(&mut self, author: Author) {
        self.authors.push(author);
    }

    pub fn set_doi(&mut self, doi: String) {
        self.doi = Some(doi);
    }

    pub fn set_pmid(&mut self, pmid: String) {
        self.pmid = Some(pmid);
    }

    pub fn set_pmcid(&mut self, pmcid: String) {
        self.pmcid = Some(pmcid);
    }

    pub fn set_year(&mut self, year: String) {
        self.year = Some(year);
    }

    pub fn set_journal(&mut self, journal: String) {
        self.journal = Some(journal);
    }

    pub fn build(self) -> Result<Reference, ReferenceBuildError> {
        let title = self.title.ok_or_else(|| ReferenceBuildError::NoTitle)?;
        let external = RefExternal {
            pmid: self.pmid,
            pmcid: self.pmcid,
            doi: self.doi,
        };

        Ok(Reference {
            title,
            authors: self.authors,
            journal: self.journal.unwrap(),
            year: self.year.unwrap(),
            external,
        })
    }

    pub fn clear(&mut self) {
        self.title = None;
        self.authors.clear();
        self.journal = None;
        self.year = None;
        self.pmid = None;
        self.pmcid = None;
        self.doi = None;
    }
}
