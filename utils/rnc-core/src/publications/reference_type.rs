use std::{
    fmt,
    str::FromStr,
};

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ConversionError {
    #[error("The prefix `{0}` is not known")]
    UnknownReferenceType(String),
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, Hash)]
pub enum ReferenceType {
    Pmid,
    Doi,
    Pmcid,
}

impl<'a> From<ReferenceType> for &'a str {
    fn from(raw: ReferenceType) -> &'a str {
        match raw {
            ReferenceType::Pmid => "pmid",
            ReferenceType::Doi => "doi",
            ReferenceType::Pmcid => "pmcid",
        }
    }
}

impl FromStr for ReferenceType {
    type Err = ConversionError;

    fn from_str(raw: &str) -> Result<ReferenceType, Self::Err> {
        match raw {
            "pmid" => Ok(ReferenceType::Pmid),
            "pmcid" => Ok(ReferenceType::Pmcid),
            "doi" => Ok(ReferenceType::Doi),
            _ => Err(ConversionError::UnknownReferenceType(raw.to_string())),
        }
    }
}

impl fmt::Display for ReferenceType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Pmid => write!(f, "pmid"),
            Self::Doi => write!(f, "doi"),
            Self::Pmcid => write!(f, "pmcid"),
        }
    }
}
