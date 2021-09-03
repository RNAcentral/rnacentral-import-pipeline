use serde::{
    de::DeserializeOwned,
    Deserialize,
    Serialize,
};

use std::{
    cmp::Ordering::{
        Equal,
        Greater,
        Less,
    },
    path::Path,
};

use anyhow::{
    anyhow,
    Result,
};
use itertools::Itertools;

use crate::psql::PsqlJsonIterator;
use rnc_utils;

pub trait HasIndex {
    fn index(&self) -> usize;
}

pub enum Criteria {
    AnyNumber,
    AtleastOne,
    ExactlyOne,
    ZeroOrOne,
}

#[derive(Debug, Serialize, Deserialize)]
pub enum Grouped<T> {
    Multiple {
        id: usize,
        data: Vec<T>,
    },
    Optional {
        id: usize,
        data: Option<T>,
    },
    Required {
        id: usize,
        data: T,
    },
}

impl<T> Grouped<T>
where
    T: Serialize + DeserializeOwned,
{
    fn empty(criteria: &Criteria, id: usize) -> Result<Grouped<T>> {
        match criteria {
            Criteria::AnyNumber => Ok(Self::Multiple {
                id,
                data: Vec::new(),
            }),
            Criteria::AtleastOne => Err(anyhow!("Missing data for id: {}", id)),
            Criteria::ExactlyOne => Err(anyhow!("Missing data for id: {}", id)),
            Criteria::ZeroOrOne => Ok(Self::Optional {
                id,
                data: None,
            }),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            Self::Multiple { id: _id, data } => data.len() == 0,
            Self::Optional { id: _id, data } => data.is_none(),
            Self::Required { id: _id, data: _data } => false,
        }
    }

    fn build(criteria: &Criteria, id: usize, items: Vec<T>) -> Result<Grouped<T>> {
        match (criteria, items.len()) {
            (Criteria::AnyNumber, 0) => Ok(Self::Multiple {
                id,
                data: Vec::new(),
            }),
            (Criteria::AnyNumber, _) => Ok(Self::Multiple {
                id,
                data: items,
            }),
            (Criteria::AtleastOne, 0) => Err(anyhow!("Missing at least one item for {}", id)),
            (Criteria::AtleastOne, _) => Ok(Self::Multiple {
                id,
                data: items,
            }),
            (Criteria::ExactlyOne, 0) => Err(anyhow!("Missing required items for {}", id)),
            (Criteria::ExactlyOne, 1) => Ok(Self::Required {
                id,
                data: items.into_iter().next().unwrap(),
            }),
            (Criteria::ExactlyOne, l) => {
                Err(anyhow!("Too many items ({}) for id: {}, expected 1", l, id))
            },
            (Criteria::ZeroOrOne, 0) => Ok(Self::Optional {
                id,
                data: None,
            }),
            (Criteria::ZeroOrOne, 1) => Ok(Self::Optional {
                id,
                data: items.into_iter().next(),
            }),
            (Criteria::ZeroOrOne, l) => {
                Err(anyhow!("Too many items ({}) for {}, expected 0 or 1", l, id))
            },
        }
    }
}

pub fn group<T>(
    criteria: Criteria,
    path: &Path,
    min: usize,
    max: usize,
    output: &Path,
) -> Result<()>
where
    T: DeserializeOwned + HasIndex + Serialize,
{
    let reader = rnc_utils::buf_reader(path)?;
    let mut writer = rnc_utils::buf_writer(output)?;

    let data = PsqlJsonIterator::from_read(reader);
    let data = data.group_by(T::index);

    let mut expected = min;
    for (id, entries) in &data {
        match (&expected).cmp(&id) {
            Less => {
                while expected < id {
                    let empty: Grouped<T> = Grouped::empty(&criteria, expected)?;
                    serde_json::to_writer(&mut writer, &empty)?;
                    writeln!(&mut writer)?;
                    expected += 1;
                }
                assert!(expected == id);
                let group = Grouped::build(&criteria, expected, entries.collect())?;
                serde_json::to_writer(&mut writer, &group)?;
                writeln!(&mut writer)?;
                expected += 1;
            },
            Equal => {
                let data: Vec<T> = entries.collect();
                let grouped = Grouped::build(&criteria, id, data)?;
                serde_json::to_writer(&mut writer, &grouped)?;
                writeln!(&mut writer)?;
                expected += 1;
            },
            Greater => {
                return Err(anyhow!(
                    "Somehow got too small index {}, expected, {}",
                    &id,
                    &expected
                ));
            },
        }

        if expected != id + 1 {
            return Err(anyhow!("Current out of sync with expected {} {}", expected, id));
        }
    }

    match expected.cmp(&max) {
        Less => {
            while expected < max {
                let empty: Grouped<T> = Grouped::empty(&criteria, expected)?;
                serde_json::to_writer(&mut writer, &empty)?;
                writeln!(&mut writer)?;
                expected += 1;
            }
        },
        Equal => (),
        Greater => {
            return Err(anyhow!("Got more items {} than expected {}", &expected, &max));
        },
    }

    assert!(expected == max, "Did not reach {}, at {}", &max, &expected);

    Ok(())
}
