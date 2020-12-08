use std::str::FromStr;

use thiserror::Error;

use crate::publications::reference_type;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct ExternalReference(reference_type::ReferenceType, String);

#[derive(Error, Debug)]
pub enum ConversionError {
    #[error("No prefix for the reference id: `{0}`")]
    MissingPrefix(String),

    #[error("Unknown type of reference: {0}")]
    RefTypeError(#[from] reference_type::ConversionError),

    #[error("Format of reference `{0}` is invalid")]
    InvalidFormat(String),
}

impl ExternalReference {
    pub fn new(ref_type: reference_type::ReferenceType, ref_id: String) -> Self {
        Self(ref_type, ref_id)
    }

    pub fn pmid(ref_id: String) -> Self {
        Self(reference_type::ReferenceType::Pmid, ref_id)
    }

    pub fn pmcid(ref_id: String) -> Self {
        Self(reference_type::ReferenceType::Pmcid, ref_id)
    }

    pub fn doi(ref_id: String) -> Self {
        Self(reference_type::ReferenceType::Doi, ref_id)
    }

    pub fn ref_type(&self) -> reference_type::ReferenceType {
        self.0
    }

    pub fn ref_id(&self) -> &str {
        &self.1
    }

    pub fn to_string(&self) -> String {
        format!("{}:{}", self.0.to_string(), self.0)
    }
}

impl FromStr for ExternalReference {
    type Err = ConversionError;

    fn from_str(raw: &str) -> Result<ExternalReference, Self::Err> {
        let parts: Vec<&str> = raw.split(":").collect();
        if parts.len() == 1 {
            return Err(Self::Err::MissingPrefix(raw.to_string()));
        }

        if parts.len() > 2 {
            return Err(Self::Err::InvalidFormat(raw.to_string()));
        }

        let ref_type: reference_type::ReferenceType = parts[0].parse()?;
        Ok(ExternalReference(ref_type, parts[1].to_string()))
    }
}

impl From<ExternalReference> for String {
    fn from(raw: ExternalReference) -> String {
        format!("{}:{}", raw.0, raw.1)
    }
}
