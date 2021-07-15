use std::{
    fs::File,
    io,
    io::{
        prelude::*,
        BufReader,
        BufWriter,
    },
    path::Path,
};

pub fn buf_reader(filename: &Path) -> io::Result<Box<dyn BufRead>> {
    if filename == Path::new("-") {
        let stdin = io::stdin();
        stdin.lock();
        return Ok(Box::new(BufReader::new(stdin)));
    }
    let file = File::open(filename)?;
    Ok(Box::new(BufReader::new(file)))
}

pub fn buf_writer(filename: &Path) -> io::Result<Box<dyn Write>> {
    if filename == Path::new("-") {
        return Ok(Box::new(BufWriter::new(io::stdout())));
    }
    let file = File::create(filename)?;
    Ok(Box::new(BufWriter::new(file)))
}

pub fn writer(filename: &Path) -> io::Result<Box<dyn Write>> {
    if filename == Path::new("-") {
        return Ok(Box::new(io::stdout()));
    }
    let file = File::create(filename)?;
    Ok(Box::new(file))
}

pub fn reader(filename: &Path) -> io::Result<Box<dyn Read>> {
    if filename == Path::new("-") {
        let stdin = io::stdin();
        stdin.lock();
        return Ok(Box::new(stdin));
    }
    let file = File::open(filename)?;
    Ok(Box::new(file))
}
