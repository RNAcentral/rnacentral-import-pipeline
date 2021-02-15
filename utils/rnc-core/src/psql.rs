use std::{
    fs::File,
    io::{
        BufRead,
        BufReader,
        Read,
    },
    path::Path,
};

use std::marker::PhantomData;

use anyhow::{
    Context,
    Result,
};

use serde::de::DeserializeOwned;

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
    pub fn from_path(path: &Path) -> Result<Self> {
        let file = File::open(path).with_context(|| format!("Could not open file {:?}", &path))?;
        Ok(JsonlIterator::from_read(file))
    }
}

impl<R: Read, T: DeserializeOwned> Iterator for JsonlIterator<R, T> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        match self.reader.read_line(&mut self.buf).unwrap() {
            0 => None,
            _ => {
                let value: T = serde_json::from_str(&self.buf)
                    .or(serde_json::from_str(&self.buf.replace("////", "//")))
                    .with_context(|| format!("Could not parse line {}", &self.buf))
                    .unwrap();
                self.buf.clear();
                Some(value)
            },
        }
    }
}
