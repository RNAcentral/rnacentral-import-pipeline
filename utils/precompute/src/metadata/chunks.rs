use serde::Deserialize;

use std::{
    fs::File,
    io::BufReader,
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

use csv::Reader;

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
                end: raw.stop,
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
