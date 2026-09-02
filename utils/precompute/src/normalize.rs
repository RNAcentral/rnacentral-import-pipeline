use std::{
    cell::Cell,
    path::Path,
};

use serde::{
    Deserialize,
    Serialize,
};

use anyhow::Result;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use crate::{
    accessions::{
        Accession,
        RawAccessionEntry,
    },
    metadata::{
        coordinate::Coordinate,
        merged::Metadata,
        orf::OrfInfo,
        previous::Previous,
        r2dt_hit::R2dtHit,
        rfam_hit::RfamHit,
    },
};

use rnc_core::{
    grouper::Grouped,
    psql::PsqlJsonIterator,
};

#[derive(Debug, Deserialize, Serialize)]
pub struct Normalized {
    upi: String,
    taxid: usize,
    length: usize,
    last_release: usize,
    coordinates: Vec<Coordinate>,
    accessions: Vec<Accession>,
    deleted: bool,
    previous: Option<Previous>,
    rfam_hits: Vec<RfamHit>,
    r2dt_hits: Vec<R2dtHit>,
    orf_info: Option<OrfInfo>,
    possible_orf: Option<bool>,
    possible_orf_stopfree: Option<bool>,
    possible_orf_tcode: Option<bool>,
}

impl Normalized {
    fn new(raw_accessions: Vec<RawAccessionEntry>, metadata: Metadata) -> Result<Self> {
        assert!(raw_accessions.len() != 0, "Must given accessions to normalize");
        let accessions: Vec<Accession> = raw_accessions.into_iter().map(Accession::from).collect();
        let last_release = accessions.iter().map(|a| a.last_release).max().unwrap();
        let deleted = accessions.iter().all(|a| !a.is_active);

        return Ok(Self {
            upi: metadata.upi,
            taxid: metadata.taxid,
            length: metadata.length,
            last_release,
            coordinates: metadata.coordinates,
            accessions,
            deleted,
            previous: metadata.previous,
            rfam_hits: metadata.rfam_hits,
            r2dt_hits: metadata.r2dt_hits.into_iter().collect(),
            orf_info: metadata.orf_info,
            possible_orf: metadata.possible_orf,
            possible_orf_stopfree: metadata.possible_orf_stopfree,
            possible_orf_tcode: metadata.possible_orf_tcode,
        });
    }
}

pub fn write(accession_file: &Path, metadata_file: &Path, output: &Path) -> Result<()> {
    let accessions = PsqlJsonIterator::from_path(accession_file)?;
    let accessions = accessions.map(|group: Grouped<RawAccessionEntry>| match group {
        Grouped::Multiple {
            id,
            data,
        } => (id, data),
        _ => panic!("Illegal data format for accessions file {:?}", &group),
    });

    // A pair with no accessions left cannot be described, so skip it instead of
    // tripping the assert in Normalized::new.
    let skipped = Cell::new(0usize);
    let accessions = accessions.filter(|(_id, data)| {
        if data.is_empty() {
            skipped.set(skipped.get() + 1);
            return false;
        }
        true
    });
    let accessions = accessions.into_iter().assume_sorted_by_key();

    let metadata = PsqlJsonIterator::from_path(metadata_file)?;
    let metadata = metadata.map(|m: Metadata| (m.id, m));
    let metadata = metadata.into_iter().assume_sorted_by_key();
    let partial = accessions.join(metadata);

    let mut output = rnc_utils::buf_writer(output)?;
    for (_id, data) in partial {
        let (accessions, metadata) = data;
        let norm = Normalized::new(accessions, metadata)?;
        serde_json::to_writer(&mut output, &norm)?;
        writeln!(&mut output)?;
    }

    if skipped.get() > 0 {
        log::warn!("Skipped {} urs_taxid with no accessions", skipped.get());
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs,
        path::PathBuf,
        sync::atomic::{
            AtomicUsize,
            Ordering as AtomicOrdering,
        },
    };

    fn temp_path(name: &str) -> PathBuf {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, AtomicOrdering::Relaxed);
        std::env::temp_dir().join(format!(
            "precompute-normalize-{}-{}-{}",
            std::process::id(),
            n,
            name
        ))
    }

    fn accession_group(id: usize, count: usize) -> String {
        let entries: Vec<String> = (0..count)
            .map(|n| {
                format!(
                    r#"{{"id":{id},"urs_id":{id},"urs_taxid":"URS00000000{id:02}_9606","accession":"A{id}{n}","last_release":1,"is_active":true,"description":"An RNA","gene":null,"optional_id":null,"database":"ENA","species":null,"common_name":null,"feature_name":null,"ncrna_class":null,"locus_tag":null,"organelle":null,"lineage":null,"all_species":[],"all_common_names":[],"so_rna_type":null}}"#,
                    id = id,
                    n = n
                )
            })
            .collect();
        format!(r#"{{"Multiple":{{"id":{},"data":[{}]}}}}"#, id, entries.join(","))
    }

    fn metadata(id: usize) -> String {
        format!(
            r#"{{"id":{id},"urs_id":{id},"urs_taxid":"URS00000000{id:02}_9606","upi":"URS00000000{id:02}","taxid":9606,"length":100,"coordinates":[],"previous":null,"rfam_hits":[],"r2dt_hits":null,"orf_info":null,"possible_orf":null,"possible_orf_stopfree":null,"possible_orf_tcode":null}}"#,
            id = id
        )
    }

    /// An empty accession group is dropped instead of tripping the
    /// "Must given accessions to normalize" assert; the pairs around it still
    /// come through.
    #[test]
    fn skips_urs_taxid_without_accessions() -> Result<()> {
        let accessions = temp_path("accessions.json");
        let meta = temp_path("metadata.json");
        let output = temp_path("merged.json");

        fs::write(
            &accessions,
            format!(
                "{}\n{}\n{}\n",
                accession_group(1, 1),
                accession_group(2, 0),
                accession_group(3, 2)
            ),
        )?;
        fs::write(&meta, format!("{}\n{}\n{}\n", metadata(1), metadata(2), metadata(3)))?;

        write(&accessions, &meta, &output)?;

        let written = fs::read_to_string(&output)?;
        let upis: Vec<String> = written
            .lines()
            .map(|l| serde_json::from_str::<serde_json::Value>(l).unwrap()["upi"].to_string())
            .collect();
        assert_eq!(upis, vec!["\"URS0000000001\"", "\"URS0000000003\""]);

        fs::remove_file(&accessions)?;
        fs::remove_file(&meta)?;
        fs::remove_file(&output)?;
        Ok(())
    }
}
