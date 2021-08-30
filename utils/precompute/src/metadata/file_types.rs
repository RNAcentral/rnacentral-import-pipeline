use std::{
    str::FromStr,
    string::ToString,
};
use strum_macros;

#[derive(strum_macros::ToString, strum_macros::EnumString, Debug)]
pub enum FileType {
    #[strum(serialize = "basic")]
    Basic,
    #[strum(serialize = "coordinates")]
    Coordinates,
    #[strum(serialize = "orfs")]
    Orfs,
    #[strum(serialize = "previous")]
    Previous,
    #[strum(serialize = "r2dt-hits")]
    R2dtHits,
    #[strum(serialize = "rfam-hits")]
    RfamHits,
}
