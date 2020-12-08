use std::{
    str,
    str::FromStr,
};

use thiserror::Error;

use crate::urs::Urs;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Could not parse urs: {0}")]
    CannotParseUrs(String),

    #[error("Could not parse taxid: {0}")]
    CannotParseTaxid(String),
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct UrsTaxid(u64, u64);

impl UrsTaxid {
    pub fn new(urs: u64, taxid: u64) -> Self {
        UrsTaxid(urs, taxid)
    }

    pub fn to_string(&self) -> String {
        format!("URS{:010X}_{}", self.0, self.1)
    }

    pub fn urs(&self) -> u64 {
        self.0
    }

    pub fn taxid(&self) -> u64 {
        self.1
    }

    pub fn short(&self) -> String {
        let urs: Urs = self.urs().into();
        format!("{}_{}", urs.short_urs(), self.taxid())
    }
}

impl FromStr for UrsTaxid {
    type Err = Error;

    fn from_str(raw: &str) -> Result<Self, Self::Err> {
        let (raw_urs, raw_taxid) = raw.split_at(14);
        let urs = u64::from_str_radix(&raw_urs[3..13], 16)
            .map_err(|_| Error::CannotParseUrs(raw.to_string()))?;

        let taxid =
            raw_taxid.parse::<u64>().map_err(|_| Error::CannotParseTaxid(raw_taxid.to_string()))?;
        Ok(Self(urs, taxid))
    }
}

impl From<&UrsTaxid> for String {
    fn from(urs: &UrsTaxid) -> String {
        format!("URS{:010X}_{}", urs.0, urs.1)
    }
}

impl From<&UrsTaxid> for Urs {
    fn from(urs: &UrsTaxid) -> Urs {
        Urs::from(urs.0)
    }
}

impl From<UrsTaxid> for Urs {
    fn from(urs: UrsTaxid) -> Urs {
        Urs::from(urs.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;

    #[test]
    fn can_convert_string_to_urs_taxid() -> Result<(), Box<dyn Error>> {
        assert_eq!("URS0000000009_1".parse::<UrsTaxid>()?, UrsTaxid(9, 1));
        assert_eq!("URS0000C0472E_12445".parse::<UrsTaxid>()?, UrsTaxid(12601134, 12445));
        assert_eq!("URS0000000001_562".parse::<UrsTaxid>()?, UrsTaxid(1, 562));
        assert_eq!("URS00008B8A75_9606".parse::<UrsTaxid>()?, UrsTaxid(9144949, 9606));
        assert_eq!("URS00001EE391_1250050".parse::<UrsTaxid>()?, UrsTaxid(2024337, 1250050));
        assert_eq!("URS00008C3642_9606".parse::<UrsTaxid>()?, UrsTaxid(9188930, 9606));
        Ok(())
    }
}
