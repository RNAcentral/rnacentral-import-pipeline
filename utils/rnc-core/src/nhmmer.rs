use regex::Regex;

pub fn valid_sequence(sequence: &str) -> bool {
    lazy_static! {
        static ref RE: Regex = Regex::new(r"^[ACGTN]+$").unwrap();
    }
    RE.is_match(sequence)
}
