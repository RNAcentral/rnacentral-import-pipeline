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
            Self::Multiple {
                id: _id,
                data,
            } => data.len() == 0,
            Self::Optional {
                id: _id,
                data,
            } => data.is_none(),
            Self::Required {
                id: _id,
                data: _data,
            } => false,
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

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::{
        json,
        Value,
    };
    use std::{
        fs,
        io::Write,
    };

    #[derive(Debug, Serialize, Deserialize)]
    struct Item {
        id: usize,
        name: String,
    }

    impl HasIndex for Item {
        fn index(&self) -> usize {
            self.id
        }
    }

    fn run(criteria: Criteria, ids: &[usize], min: usize, max: usize) -> Result<Vec<Value>> {
        let dir = tempfile::tempdir()?;
        let input = dir.path().join("in.json");
        let output = dir.path().join("out.json");
        let mut file = fs::File::create(&input)?;
        for id in ids {
            writeln!(file, "{}", json!({"id": id, "name": format!("n{}", id)}))?;
        }
        drop(file);

        group::<Item>(criteria, &input, min, max, &output)?;
        let written = fs::read_to_string(&output)?;
        Ok(written.lines().map(|l| serde_json::from_str(l).unwrap()).collect())
    }

    fn item(id: usize) -> Value {
        json!({"id": id, "name": format!("n{}", id)})
    }

    #[test]
    fn any_number_fills_the_gaps_with_empty_groups() {
        let out = run(Criteria::AnyNumber, &[1, 3], 1, 5).unwrap();
        assert_eq!(
            out,
            vec![
                json!({"Multiple": {"id": 1, "data": [item(1)]}}),
                json!({"Multiple": {"id": 2, "data": []}}),
                json!({"Multiple": {"id": 3, "data": [item(3)]}}),
                json!({"Multiple": {"id": 4, "data": []}}),
            ]
        );
    }

    #[test]
    fn rows_with_one_id_are_collected_into_one_group() {
        let out = run(Criteria::AtleastOne, &[1, 1, 2], 1, 3).unwrap();
        assert_eq!(out[0], json!({"Multiple": {"id": 1, "data": [item(1), item(1)]}}));
        assert_eq!(out.len(), 2);
    }

    #[test]
    fn zero_or_one_leaves_a_missing_id_as_none() {
        let out = run(Criteria::ZeroOrOne, &[2], 1, 3).unwrap();
        assert_eq!(
            out,
            vec![
                json!({"Optional": {"id": 1, "data": null}}),
                json!({"Optional": {"id": 2, "data": item(2)}}),
            ]
        );
    }

    #[test]
    fn exactly_one_wraps_the_single_row_as_required() {
        let out = run(Criteria::ExactlyOne, &[1, 2], 1, 3).unwrap();
        assert_eq!(out[1], json!({"Required": {"id": 2, "data": item(2)}}));
    }

    #[test]
    fn exactly_one_fails_for_a_missing_id() {
        let err = run(Criteria::ExactlyOne, &[1, 3], 1, 4).unwrap_err();
        assert!(err.to_string().contains("Missing data for id: 2"), "{}", err);
    }

    #[test]
    fn exactly_one_fails_for_a_repeated_id() {
        let err = run(Criteria::ExactlyOne, &[1, 1], 1, 2).unwrap_err();
        assert!(err.to_string().contains("Too many items"), "{}", err);
    }

    #[test]
    fn zero_or_one_fails_for_a_repeated_id() {
        let err = run(Criteria::ZeroOrOne, &[1, 1], 1, 2).unwrap_err();
        assert!(err.to_string().contains("expected 0 or 1"), "{}", err);
    }

    #[test]
    fn fails_when_the_data_runs_past_the_upper_bound() {
        let err = run(Criteria::AnyNumber, &[1, 2, 3], 1, 2).unwrap_err();
        assert!(err.to_string().contains("more items"), "{}", err);
    }

    #[test]
    fn fails_when_an_id_is_below_the_lower_bound() {
        let err = run(Criteria::AnyNumber, &[1], 2, 3).unwrap_err();
        assert!(err.to_string().contains("too small index"), "{}", err);
    }

    #[test]
    fn emptiness_follows_the_variant() {
        assert!(Grouped::<Item>::Multiple {
            id: 1,
            data: vec![]
        }
        .is_empty());
        assert!(Grouped::<Item>::Optional {
            id: 1,
            data: None
        }
        .is_empty());
        let required = Grouped::Required {
            id: 1,
            data: Item {
                id: 1,
                name: "a".into(),
            },
        };
        assert!(!required.is_empty());
    }
}
