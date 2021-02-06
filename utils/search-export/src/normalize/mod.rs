use std::path::Path;

use anyhow::Result;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use itertools::Itertools;

pub mod ds;
pub mod utils;

use crate::normalize::ds::so_tree;
use rnc_core::psql::JsonlIterator;

pub fn write_merge(
    base_file: &Path,
    crs_file: &Path,
    feedback_file: &Path,
    go_annotations_file: &Path,
    interacting_proteins_file: &Path,
    interacting_rnas_file: &Path,
    precompute_file: &Path,
    qa_status_file: &Path,
    r2dt_hits_file: &Path,
    rfam_hits_file: &Path,
    so_term_tree_file: &Path,
    output_file: &Path,
) -> Result<()> {
    let so_info = so_tree::load(&so_term_tree_file)?;
    let mut writer = rnc_utils::buf_writer(output_file)?;

    let base = JsonlIterator::from_path(base_file)?;
    let base = base.map(|b: ds::basic::Basic| (b.id, b));
    let base = base.into_iter().assume_sorted_by_key();

    let crs = JsonlIterator::from_path(crs_file)?;
    let crs = crs.group_by(|c: &ds::crs::Crs| c.id);
    let crs = crs.into_iter().assume_sorted_by_key();

    let feedback = JsonlIterator::from_path(feedback_file)?;
    let feedback = feedback.group_by(|f: &ds::feedback::Feedback| f.id);
    let feedback = feedback.into_iter().assume_sorted_by_key();

    let go_annotations = JsonlIterator::from_path(go_annotations_file)?;
    let go_annotations = go_annotations.group_by(|f: &ds::go_annotation::GoAnnotation| f.id);
    let go_annotations = go_annotations.into_iter().assume_sorted_by_key();

    let interacting_proteins = JsonlIterator::from_path(interacting_proteins_file)?;
    let interacting_proteins =
        interacting_proteins.group_by(|i: &ds::interacting_protein::InteractingProtein| i.id);
    let interacting_proteins = interacting_proteins.into_iter().assume_sorted_by_key();

    let interacting_rnas = JsonlIterator::from_path(interacting_rnas_file)?;
    let interacting_rnas =
        interacting_rnas.group_by(|i: &ds::interacting_rna::InteractingRna| i.id);
    let interacting_rnas = interacting_rnas.into_iter().assume_sorted_by_key();

    let precompute = JsonlIterator::from_path(precompute_file)?;
    let precompute = precompute.map(|p: ds::precompute::Precompute| (p.id, p));
    let precompute = precompute.into_iter().assume_sorted_by_key();

    let qa_status = JsonlIterator::from_path(qa_status_file)?;
    let qa_status = qa_status.map(|q: ds::qa_status::QaStatus| (q.id, q));
    let qa_status = qa_status.into_iter().assume_sorted_by_key();

    let r2dt_hits = JsonlIterator::from_path(r2dt_hits_file)?;
    let r2dt_hits = r2dt_hits.map(|h: ds::r2dt::R2dt| (h.id, h));
    let r2dt_hits = r2dt_hits.into_iter().assume_sorted_by_key();

    let rfam_hits = JsonlIterator::from_path(rfam_hits_file)?;
    let rfam_hits = rfam_hits.group_by(|i: &ds::rfam_hit::RfamHit| i.id);
    let rfam_hits = rfam_hits.into_iter().assume_sorted_by_key();

    let merged = base
        .join(precompute)
        .join(qa_status)
        .left_join(crs)
        .left_join(feedback)
        .left_join(go_annotations)
        .left_join(interacting_proteins)
        .left_join(interacting_rnas)
        .left_join(r2dt_hits)
        .left_join(rfam_hits);

    for (id, entry) in merged {
        let (
            (
                (
                    (
                        (((((basic, precompute), qa_status), crs), feedback), go_annotations),
                        interacting_proteins,
                    ),
                    interacting_rnas,
                ),
                r2dt_hits,
            ),
            rfam_hits,
        ) = entry;

        let so_rna_type_tree = so_info[precompute.so_rna_type()].clone();

        let raw = ds::entry::Raw {
            id,
            base: basic,
            precompute,
            qa_status,
            crs: crs.into_iter().map(|g| g.into_iter()).flatten().collect(),
            feedback: feedback.into_iter().map(|g| g.into_iter()).flatten().collect(),
            go_annotations: go_annotations.into_iter().map(|g| g.into_iter()).flatten().collect(),
            interacting_proteins: interacting_proteins
                .into_iter()
                .map(|g| g.into_iter())
                .flatten()
                .collect(),
            interacting_rnas: interacting_rnas
                .into_iter()
                .map(|g| g.into_iter())
                .flatten()
                .collect(),
            r2dt: r2dt_hits,
            rfam_hits: rfam_hits.into_iter().map(|g| g.into_iter()).flatten().collect(),
            so_tree: so_rna_type_tree,
        };

        serde_json::to_writer(&mut writer, &raw)?;
        writeln!(&mut writer)?;
    }

    Ok(())
}

pub fn write(accession_file: &Path, metadata_file: &Path, output_file: &Path) -> Result<()> {
    let mut writer = rnc_utils::buf_writer(output_file)?;

    let metadata = JsonlIterator::from_path(metadata_file)?;
    let metadata = metadata.map(|r: ds::entry::Raw| (r.id, r));
    let metadata = metadata.into_iter().assume_sorted_by_key();

    let accessions = JsonlIterator::from_path(accession_file)?;
    let accessions = accessions.group_by(|a: &ds::accession::RawAccession| a.id);
    let accessions = accessions.into_iter().assume_sorted_by_key();

    let merged = accessions.join(metadata);

    for (_id, entry) in merged {
        let (accessions, metadata) = entry;
        let normalized = ds::entry::Normalized::new(metadata, accessions.collect())?;
        serde_json::to_writer(&mut writer, &normalized)?;
        writeln!(&mut writer)?;
    }

    Ok(())
}
