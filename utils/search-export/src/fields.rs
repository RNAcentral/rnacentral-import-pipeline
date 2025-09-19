use std::fmt::Display;

pub trait EntryType {
    type Fields: strum::IntoEnumIterator + Display;
    type HierarchicalField: strum::IntoEnumIterator + Display;

    fn entry_type() -> &'static str;
}

pub struct GeneEntry;
pub struct SequenceEntry;

impl EntryType for GeneEntry {
    type Fields = GeneFields;
    type HierarchicalField = SoRnaTreeField;

    fn entry_type() -> &'static str {
        "Gene"
    }
}

impl EntryType for SequenceEntry {
    type Fields = SequenceFields;
    type HierarchicalField = SoRnaTreeField;

    fn entry_type() -> &'static str {
        "Sequence"
    }
}

#[derive(strum_macros::EnumIter, Debug, PartialEq, strum_macros::Display)]
#[strum(serialize_all = "snake_case")]
pub enum SequenceFields {
    ShortUrs,
    Length,
    Organelle,
    ExpertDb,
    CommonName,
    Function,
    Gene,
    GeneSynonym,
    #[strum(serialize = "rna_type")]
    InsdcRNAType,
    Product,
    HasGenomicCoordinates,
    Md5,
    Author,
    Journal,
    InsdcSubmission,
    PubTitle,
    PubId,
    PopularSpecies,
    Boost,
    LocusTag,
    StandardName,
    RfamFamilyName,
    RfamId,
    RfamClan,
    QcWarning,
    QcWarningFound,
    TaxString,
    InvovledIn,
    PartOf,
    Enables,
    ContributesTo,
    ColocalizesWith,
    HasGoAnnotations,
    GoAnnotationSource,
    HasInteractingProteins,
    InteractingProtein,
    HasInteractingRnas,
    HasConservedStructure,
    ConservedStructure,
    OverlapsWith,
    NoOverlapsWith,
    HasSecondaryStructure,
    SecondaryStructureModel,
    SecondaryStructureSource,
    PdbidEntityid,
    Disease,
    Url,
    OrfSource,
    HasLitScan,
    HasLitsumm,
    HasEditingEvent,
    EditChromosome,
    EditLocations,
    EditRepeatType,
    SoRnaTypeName,
    SoRnaType,
}

#[derive(strum_macros::EnumIter, Debug, PartialEq, strum_macros::Display)]
#[strum(serialize_all = "snake_case")]
pub enum GeneFields {
    Length,
    ExpertDb,
    Gene,
    PubliGeneName,
    GeneSynonym,
    #[strum(serialize = "rna_type")]
    InsdcRNAType,
    HasGoAnnotations,
    HasGenomicCoordinates,
    Boost,
    StandardName,
    QcWarningFound,
    GeneMember,
    HasSecondaryStructure,
    HasLitScan,
    HasLitsumm,
    HasEditingEvent,
    SoRnaTypeName,
}

#[derive(strum_macros::EnumIter, Debug, PartialEq, strum_macros::Display)]
#[strum(serialize_all = "snake_case")]
pub enum SoRnaTreeField {
    SoRnaType,
}
