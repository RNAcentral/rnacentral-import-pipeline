use std::{
    fs::File,
    io::{
        BufWriter,
        Write,
    },
    path::Path,
};

use anyhow::Result;

use sorted_iter::{
    assume::*,
    SortedPairIterator,
};

use itertools::Itertools;

use rnc_core::psql::JsonlIterator;

use crate::{
    genes::{
        gene::Entry,
        info::GeneInfo,
        member::GeneMembers,
    },
    sequences::entry::Normalized,
};

pub fn write_gene_members(normalized_file: &Path, output: &Path) -> Result<()> {
    let mut writer = rnc_utils::buf_writer(output)?;

    let normalized = JsonlIterator::from_path(normalized_file)?;
    let mut selected: Vec<Normalized> =
        normalized.filter(|n: &Normalized| n.is_locus_member()).collect();

    // Dangerous - this may use up too much memory one day, hopefully this stays efficent
    // as we don't have too many members. We can always split by assembly if need be
    // as this will ensure we only ever have very few sequences per chunk.
    selected.sort_by_key(|n| n.locus_id());
    let grouped = selected.into_iter().group_by(|e| e.locus_id());

    for (_id, entries) in &grouped {
        let member = GeneMembers::new(entries.collect());
        serde_json::to_writer(&mut writer, &member)?;
        writeln!(&mut writer)?;
    }

    Ok(())
}

pub fn write_gene_info(
    gene_file: &Path,
    member_file: &Path,
    xml_output: &Path,
    count_output: &Path,
) -> Result<()> {
    let genes = JsonlIterator::from_path(gene_file)?;
    let genes = genes.map(|g: GeneInfo| (g.id(), g));
    let genes = genes.into_iter().assume_sorted_by_key();

    let members = JsonlIterator::from_path(member_file)?;
    let members = members.map(|m: GeneMembers| (m.locus_id(), m));
    let members = members.into_iter().assume_sorted_by_key();

    let merged = genes.join(members);

    let mut xml_writer = BufWriter::new(File::create(&xml_output)?);
    let mut count = 0;
    for (_id, entry) in merged {
        let (info, members) = entry;
        let gene = Entry::new(info, members);
        quick_xml::se::to_writer(&mut xml_writer, &gene)?;
        writeln!(&mut xml_writer)?;
        count += 1;
    }
    xml_writer.flush()?;

    let mut count_writer = BufWriter::new(File::create(&count_output)?);
    write!(&mut count_writer, "{}", &count)?;
    count_writer.flush()?;

    Ok(())
}
