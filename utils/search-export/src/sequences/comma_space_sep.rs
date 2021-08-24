#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default)]
pub struct CommaSpaceSeparator;

impl Separator for CommaSeparator {
    #[inline]
    fn separator() -> &'static str {
        ", "
    }
}
