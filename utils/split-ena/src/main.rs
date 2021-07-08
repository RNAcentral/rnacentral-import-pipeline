use std::{
    error::Error,
    fs::{
        remove_file,
        File,
    },
    io::{
        prelude::*,
        BufReader,
        BufWriter,
    },
    path::{
        Path,
        PathBuf,
    },
};

use structopt::StructOpt;

#[derive(Debug, StructOpt)]
#[structopt(rename_all = "kebab-case")]
struct Opt {
    /// This sets the maximum number of sequences allowed in one chunk.
    #[structopt(short, long, default_value = "10000")]
    max_sequences: u64,

    /// This flag will cause this program to delete the original file once it is split
    #[structopt(short, long)]
    remove_file: bool,

    /// The name of the file to read from, using '-' means stdin.
    #[structopt(parse(from_os_str))]
    raw: PathBuf,

    /// The name of the output directory to write chunks into
    #[structopt(parse(from_os_str))]
    output_directory: PathBuf,
}

fn chunk_writer(
    directory: &Path,
    basename: &str,
    index: u64,
) -> Result<BufWriter<File>, Box<dyn Error>> {
    let mut filename = String::from(basename);
    filename.push_str("-chunk");
    filename.push_str(&index.to_string());
    let mut final_name = PathBuf::from(directory);
    final_name.push(filename);
    final_name.set_extension("ncr");

    let file = File::create(final_name)?;
    return Ok(BufWriter::new(file));
}

fn main() -> Result<(), Box<dyn Error>> {
    let opt = Opt::from_args();
    let basename = opt.raw.file_stem().unwrap().to_str().unwrap();
    let file = File::open(&opt.raw)?;
    let mut reader = BufReader::new(file);

    let mut buf = String::new();
    let mut current_count = 0;
    let mut current_chunk = 0;
    let mut writer = chunk_writer(&opt.output_directory, basename, current_chunk)?;
    loop {
        match reader.read_line(&mut buf)? {
            0 => break,
            _ => {
                writer.write_all(buf.as_bytes())?;
                if buf.starts_with("//") {
                    current_count += 1;
                    if current_count == opt.max_sequences {
                        current_count = 0;
                        current_chunk += 1;
                        writer.flush()?;
                        writer = chunk_writer(&opt.output_directory, basename, current_chunk)?;
                    }
                };
                buf.clear();
            },
        }
    }

    writer.flush()?;

    if opt.remove_file {
        remove_file(&opt.raw)?;
    }

    Ok(())
}
