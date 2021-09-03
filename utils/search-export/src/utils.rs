pub fn set_or_check<T>(maybe: &mut Option<T>, value: T)
where
    T: std::fmt::Debug + PartialEq
{
    let previous = maybe.replace(value);
    if previous.is_some() {
        assert!(&previous == maybe, "Expected {:?} and {:?} to be equal", &maybe, &previous);
    }
}
