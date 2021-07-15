use std::{
    fs::File,
    io::{
        self,
        BufRead,
        BufReader,
        Read,
    },
    marker::PhantomData,
    path::{
        Path,
        PathBuf,
    },
};

use thiserror::Error;

use serde::de::DeserializeOwned;

#[derive(Debug, Error)]
pub enum Error {
    #[error("Could not find the file {0:?}")]
    MissingFile(PathBuf),

    #[error("Could not parse line: {line}, source {source}")]
    BadlyFormattedLine {
        line: String,
        source: serde_json::Error,
    },

    #[error(transparent)]
    IoError(#[from] io::Error),
}

pub struct JsonlIterator<R: Read, T: DeserializeOwned> {
    reader: BufReader<R>,
    buf: String,
    phantom: PhantomData<T>,
}

impl<R: Read, T: DeserializeOwned> JsonlIterator<R, T> {
    pub fn from_read(reader: R) -> Self {
        Self {
            reader: BufReader::new(reader),
            buf: String::new(),
            phantom: PhantomData,
        }
    }
}

impl<T: DeserializeOwned> JsonlIterator<File, T> {
    pub fn from_path(path: &Path) -> Result<Self, Error> {
        let file = File::open(path)?;
        Ok(JsonlIterator::from_read(file))
    }
}

impl<R: Read, T: DeserializeOwned> Iterator for JsonlIterator<R, T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        match self.reader.read_line(&mut self.buf).unwrap() {
            0 => None,
            _ => {
                let value: T = serde_json::from_str(&self.buf.replace("\\\\", "\\"))
                    .or(serde_json::from_str(&self.buf))
                    .map_err(|source| Error::BadlyFormattedLine {
                        line: self.buf.to_string(),
                        source,
                    })
                    .unwrap();
                self.buf.clear();
                Some(value)
            },
        }
    }
}
