use thiserror::Error;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Expecting a non-empty {0} list")]
    UnexpectedEmpty(String),

    #[error("Expecting a single {0} value")]
    TooMany(String),
}

pub fn expect_single<T: Clone>(values: &[T], name: &str) -> Result<T, Error> {
    match values.len() {
        0 => Err(Error::UnexpectedEmpty(name.to_string())),
        1 => {
            let val = values.first().cloned().unwrap();
            Ok(val)
        },
        _ => Err(Error::TooMany(name.to_string())),
    }
}

pub fn maybe_single<T: Clone>(values: &[T], name: &str) -> Result<Option<T>, Error> {
    match values.len() {
        0 => Ok(None),
        1 => {
            let val = values.first().cloned().unwrap();
            Ok(Some(val))
        },
        _ => Err(Error::TooMany(name.to_string())),
    }
}
