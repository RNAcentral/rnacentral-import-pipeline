use std::{
    cmp::Ordering::{
        Equal,
        Greater,
        Less,
    },
    fmt::Debug,
    fs::File,
    io::BufReader,
    path::Path,
};

use serde::{
    de::DeserializeOwned,
    Deserialize,
    Serialize,
};

use itertools::Itertools;
use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use anyhow::{
    anyhow,
    Result,
};

/// A `(id, urs, release)` record. Used by the `max-release` command.
#[derive(Serialize, Deserialize, Debug)]
pub struct UrsEntry {
    id: usize,
    urs: String,
    release: usize,
}

/// A `(urs_taxid, release)` record used for release-based selection. Both the
/// xref and precompute inputs are produced per `urs_taxid` (already aggregated
/// to one row per pair by the SQL), so the `urs_taxid` is the natural join key.
#[derive(Serialize, Deserialize, Debug)]
pub struct ReleaseEntry {
    urs_taxid: String,
    release: usize,
}

fn entries<T: DeserializeOwned>(path: &Path) -> Result<impl Iterator<Item = T>> {
    let file = File::open(path)?;
    let buf = BufReader::new(file);
    let reader = csv::ReaderBuilder::new().has_headers(false).from_reader(buf);

    let records = reader.into_deserialize().filter_map(Result::ok);

    Ok(records)
}

/// Wrap an iterator of `(key, value)` and panic if the keys are ever
/// decreasing. The release selection merge-join relies on both inputs being
/// sorted by key; this turns a silent mis-join (dropped/mismatched records)
/// into a loud, immediate failure.
fn checked_sorted_by_key<I, K, V>(iter: I) -> impl Iterator<Item = (K, V)>
where
    I: Iterator<Item = (K, V)>,
    K: Ord + Clone + Debug,
{
    let mut previous: Option<K> = None;
    iter.map(move |(key, value)| {
        if let Some(ref prev) = previous {
            if key < *prev {
                panic!(
                    "release selection input is not sorted by key: {:?} came after {:?}",
                    key, prev
                );
            }
        }
        previous = Some(key.clone());
        (key, value)
    })
}

pub fn write_max(filename: &Path, output: &Path) -> Result<()> {
    let out = rnc_utils::buf_writer(output)?;
    let mut writer = csv::Writer::from_writer(out);

    let results = entries::<UrsEntry>(filename)?.group_by(|e: &UrsEntry| e.id);

    for (_, entries) in &results {
        let max = entries.max_by(|l, r| l.release.cmp(&r.release));
        match max {
            Some(v) => writer.serialize(v)?,
            None => (),
        }
    }
    writer.flush()?;

    Ok(())
}

/// Selects which `urs_taxid` pairs to forward to the precompute workflow, based
/// on when they were last updated.
///
/// # Description
///
/// `xref` is updated as a sequence is imported, so the `release` value there is
/// the last release in which we touched the `urs_taxid` (a deletion also bumps
/// it). `known` comes from the existing precompute table. Comparing the two per
/// `urs_taxid`:
///
/// * if the known release equals the xref release the pair is up to date and is skipped,
/// * if the known release is less than the xref release (or the pair has never been
///   precomputed) the pair is selected for recomputation,
/// * if the known release is greater than the xref release something has gone wrong and
///   we error out.
///
/// Both inputs must be CSV with columns `(urs_taxid, release)`, one row per
/// `urs_taxid`, sorted by `urs_taxid` (byte order, i.e. `COLLATE "C"` on the
/// SQL side so it matches Rust's string ordering). The selected `urs_taxid`
/// values are written, one per line, to `output`.
///
/// # Arguments
///
/// * `xrefs` - Path to the CSV of `(urs_taxid, release)` from `xref`.
/// * `known` - Path to the CSV of `(urs_taxid, release)` from the precompute table.
/// * `output` - Path to write the selected `urs_taxid` values to.
///
/// # Errors
///
/// Returns an error if a file cannot be read/written, or if any precompute
/// release is newer than the corresponding xref release.
pub fn select_new(xrefs: &Path, known: &Path, output: &Path) -> Result<()> {
    let xref_records =
        checked_sorted_by_key(entries::<ReleaseEntry>(xrefs)?.map(|e| (e.urs_taxid.clone(), e)))
            .assume_sorted_by_key();
    let known_records =
        checked_sorted_by_key(entries::<ReleaseEntry>(known)?.map(|e| (e.urs_taxid.clone(), e)))
            .assume_sorted_by_key();

    let mut writer = csv::Writer::from_writer(File::create(output)?);
    let pairs = xref_records.outer_join(known_records);
    for (_key, (xref, pre)) in pairs {
        match (xref, pre) {
            (Some(x), Some(p)) => match x.release.cmp(&p.release) {
                Less => Err(anyhow!(
                    "This should never happen, too small release for {:?} vs {:?}",
                    &x,
                    &p
                ))?,
                Equal => (),
                Greater => writer.write_record(&[x.urs_taxid])?,
            },
            (Some(x), None) => writer.write_record(&[x.urs_taxid])?,
            (None, Some(_)) => (),
            (None, None) => (),
        }
    }
    writer.flush()?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        io::{
            BufRead,
            Write,
        },
        sync::atomic::{
            AtomicUsize,
            Ordering as AtomicOrdering,
        },
    };

    /// Unique-per-run temp filename. Process id plus a counter is enough here
    /// and keeps `rand` out of the dependency tree.
    fn get_random_fname() -> String {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, AtomicOrdering::Relaxed);
        std::env::temp_dir()
            .join(format!("precompute-releases-{}-{}.csv", std::process::id(), n))
            .to_string_lossy()
            .into_owned()
    }

    /// Write `contents` to a fresh temp file and return its path.
    fn write_temp(contents: &str) -> String {
        let fname = get_random_fname();
        let mut file = File::create(&fname).unwrap();
        file.write_all(contents.as_bytes()).unwrap();
        fname
    }

    /// Read the output file into a sorted Vec of lines.
    fn read_sorted_lines(path: &Path) -> Vec<String> {
        let file = File::open(path).unwrap();
        let mut lines: Vec<String> = BufReader::new(file).lines().map(|l| l.unwrap()).collect();
        lines.sort();
        lines
    }

    #[test]
    fn test_select_new() -> Result<()> {
        // One row per urs_taxid per side, sorted by urs_taxid (byte order).
        // URS1: xref newer -> selected
        // URS2: equal -> skipped
        // URS3: xref newer -> selected
        // URS6: only in xref (novel) -> selected
        // URS7: only in known (gone from xref) -> skipped
        let xref_data = "URS1_9606,604\n\
                         URS2_9606,600\n\
                         URS3_9606,605\n\
                         URS6_9606,607\n";
        let known_data = "URS1_9606,600\n\
                          URS2_9606,600\n\
                          URS3_9606,600\n\
                          URS7_9606,600\n";

        let xref_fname = write_temp(xref_data);
        let known_fname = write_temp(known_data);
        let output_fname = get_random_fname();
        let output = Path::new(&output_fname);

        select_new(Path::new(&xref_fname), Path::new(&known_fname), output)?;

        let got = read_sorted_lines(output);
        assert_eq!(
            got,
            vec!["URS1_9606".to_string(), "URS3_9606".to_string(), "URS6_9606".to_string(),]
        );

        std::fs::remove_file(output)?;
        std::fs::remove_file(Path::new(&xref_fname))?;
        std::fs::remove_file(Path::new(&known_fname))?;

        Ok(())
    }

    #[test]
    fn test_select_new_known_newer() -> Result<()> {
        // precompute release newer than xref -> catastrophic, must error.
        let xref_data = "URS1_9606,600\n";
        let known_data = "URS1_9606,700\n";

        let xref_fname = write_temp(xref_data);
        let known_fname = write_temp(known_data);
        let output_fname = get_random_fname();
        let output = Path::new(&output_fname);

        let result = select_new(Path::new(&xref_fname), Path::new(&known_fname), output);
        assert!(result.is_err());

        let _ = std::fs::remove_file(output);
        std::fs::remove_file(Path::new(&xref_fname))?;
        std::fs::remove_file(Path::new(&known_fname))?;

        Ok(())
    }

    #[test]
    fn test_unsorted_input_panics() {
        // xref keys out of order -> the sortedness guard must fire.
        let xref_data = "URS3_9606,604\n\
                         URS1_9606,604\n";
        let known_data = "URS1_9606,600\n\
                          URS3_9606,600\n";

        let xref_fname = write_temp(xref_data);
        let known_fname = write_temp(known_data);
        let output_fname = get_random_fname();
        let output = Path::new(&output_fname);

        let result = std::panic::catch_unwind(|| {
            select_new(Path::new(&xref_fname), Path::new(&known_fname), output)
        });
        assert!(result.is_err(), "expected a panic on unsorted input");

        let _ = std::fs::remove_file(output);
        let _ = std::fs::remove_file(Path::new(&xref_fname));
        let _ = std::fs::remove_file(Path::new(&known_fname));
    }
}
