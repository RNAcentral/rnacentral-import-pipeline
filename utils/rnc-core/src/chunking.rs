use serde::Deserialize;

use std::{
    fmt,
    fs::File,
    io::{
        BufReader,
        BufWriter,
        Read,
        Write,
    },
    ops::Range,
    path::{
        Path,
        PathBuf,
    },
};

use anyhow::{
    Context,
    Result,
};

use serde::{
    de::DeserializeOwned,
    ser::Serialize,
};

use csv::Reader;

use itertools::Itertools;

use crate::psql::PsqlJsonIterator;

pub struct SingleChunk {
    endpoints: Range<usize>,
    filename: String,
}

#[derive(Debug, Deserialize)]
pub struct RawChunk {
    filename: String,
    start: usize,
    stop: usize,
}

pub struct ChunkSpec {
    ranges: Vec<SingleChunk>,
}

impl From<RawChunk> for SingleChunk {
    fn from(raw: RawChunk) -> SingleChunk {
        return Self {
            endpoints: Range {
                start: raw.start,
                end: raw.stop + 1,
            },
            filename: raw.filename,
        };
    }
}

impl SingleChunk {
    pub fn contains(&self, index: &usize) -> bool {
        self.endpoints.contains(index)
    }
}

impl ChunkSpec {
    pub fn filename_of(&self, index: usize) -> Option<PathBuf> {
        self.ranges
            .iter()
            .find(|r| r.contains(&index))
            .map(|c: &SingleChunk| PathBuf::from(&c.filename))
    }
}

pub fn load(filename: &Path) -> Result<ChunkSpec> {
    let file = File::open(filename)
        .with_context(|| format!("Could not open chunk file {:?}", &filename))?;
    let reader = BufReader::new(file);
    let mut records = Reader::from_reader(reader);
    let mut ranges = Vec::new();
    for result in records.deserialize() {
        let entry: RawChunk = result?;
        ranges.push(SingleChunk::from(entry));
    }

    Ok(ChunkSpec {
        ranges,
    })
}

pub fn write_splits<R: Read, T: DeserializeOwned + fmt::Debug + Serialize>(
    iterator: PsqlJsonIterator<R, T>,
    id_getter: fn(&T) -> usize,
    chunks: &Path,
    output: &Path,
) -> Result<()> {
    let chunks = load(chunks)?;
    let grouped = iterator.group_by(|entry: &T| {
        chunks
            .filename_of(id_getter(entry))
            .ok_or(format!("Could not find filename for {:?}", &entry))
            .unwrap()
    });

    for (filename, items) in &grouped {
        let mut path = PathBuf::from(output);
        path.push(filename);
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        for item in items {
            serde_json::to_writer(&mut writer, &item)?;
            writeln!(&mut writer)?;
        }
    }

    Ok(())
}
