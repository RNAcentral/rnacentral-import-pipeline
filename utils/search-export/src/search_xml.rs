use std::{
    collections::HashSet,
    iter::FromIterator,
    option,
    vec,
};

use log;
use serde::Serialize;
use strum::IntoEnumIterator;

use crate::{
    fields::EntryType,
    sequences::accession::CrossReference,
};

/// This trait represents the logic for converting a rust value to the XML values. This is
/// used rather an a simple Into<String>, because for some cases there is a specific
/// format.
pub trait AsSearchValue {
    fn search_value(&self) -> String;
}

impl AsSearchValue for bool {
    fn search_value(&self) -> String {
        match self {
            false => "False".to_string(),
            true => "True".to_string(),
        }
    }
}

impl AsSearchValue for String {
    fn search_value(&self) -> String {
        self.clone()
    }
}

impl AsSearchValue for &str {
    fn search_value(&self) -> String {
        self.to_string()
    }
}

impl AsSearchValue for &String {
    fn search_value(&self) -> String {
        self.to_string()
    }
}

impl AsSearchValue for usize {
    fn search_value(&self) -> String {
        self.to_string()
    }
}

enum SearchValueInner {
    Many(Vec<String>),
    Single(String),
}

pub struct SearchValue(SearchValueInner);

pub enum SearchValueIter {
    Many(vec::IntoIter<String>),
    Once(option::IntoIter<String>),
}

impl SearchValue {
    pub fn empty() -> Self {
        Self(SearchValueInner::Many(Vec::new()))
    }
}

impl<T> From<T> for SearchValue
where
    T: AsSearchValue,
{
    fn from(value: T) -> Self {
        Self(SearchValueInner::Single(value.search_value()))
    }
}

impl<T> From<Vec<T>> for SearchValue
where
    T: AsSearchValue,
{
    fn from(value: Vec<T>) -> Self {
        Self(SearchValueInner::Many(value.into_iter().map(|s| s.search_value()).collect()))
    }
}

impl<T> From<&[T]> for SearchValue
where
    T: AsSearchValue,
{
    fn from(value: &[T]) -> Self {
        // value.to_vec().into()
        Self(SearchValueInner::Many(value.iter().map(|s| s.search_value()).collect()))
    }
}

impl<T> From<HashSet<T>> for SearchValue
where
    T: AsSearchValue,
{
    fn from(value: HashSet<T>) -> Self {
        Self(SearchValueInner::Many(value.into_iter().map(|s| s.search_value()).collect()))
    }
}

impl<T> FromIterator<T> for SearchValue
where
    T: AsSearchValue,
{
    fn from_iter<A: IntoIterator<Item = T>>(iter: A) -> Self {
        let items: Vec<T> = iter.into_iter().collect();
        items.into()
    }
}

impl IntoIterator for SearchValue {
    type IntoIter = SearchValueIter;
    type Item = String;

    fn into_iter(self) -> Self::IntoIter {
        match self.0 {
            SearchValueInner::Many(vec) => SearchValueIter::Many(vec.into_iter()),
            SearchValueInner::Single(item) => SearchValueIter::Once(Some(item).into_iter()),
        }
    }
}

impl Iterator for SearchValueIter {
    type Item = String;

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Many(ref mut iter) => iter.next(),
            Self::Once(ref mut iter) => iter.next(),
        }
    }
}

pub trait SearchEntry<E>: Sized
where
    E: EntryType,
{
    /// This produces the globally unique identifier for this entry. For example, this
    /// could be a URS_taxid: `URS00000478B7_9606`.
    fn id(&self) -> &str;

    /// This produces human readable 'name'.
    fn name(&self) -> &str;

    /// A human readable description of the entry.
    fn description(&self) -> &str;

    /// Fetch the taxid for this entry. Should be something like `9606`, ie the NCBI
    /// Taxonomic id for humans.
    fn taxid(&self) -> usize;

    /// Given a field, return all values for it.
    fn field_values(&self, field: &E::Fields) -> SearchValue;

    /// Given a `Entry::HierarchicalField` this will get the value for it.
    fn tree_field_values(&self, field: &E::HierarchicalField) -> (String, Vec<String>);

    /// Fetch all cross references for this item.
    fn cross_references(&self) -> impl IntoIterator<Item = CrossReference>;

    fn search_entry(self) -> Entry {
        let mut entry = Entry::new(self.id().to_string(), self.name(), self.description());

        // All entries are active
        entry.add_field("active", true);

        // All entries have an `entry_type` value
        entry.add_field("entry_type", E::entry_type());

        // Get the values for each field
        let mut seen_any_field = false;
        for field in E::Fields::iter() {
            let name = field.to_string();
            let mut seen_field = false;
            for value in self.field_values(&field) {
                entry.add_field(&name, value);
                seen_any_field = true;
                seen_field = true;
            }

            // It should be rare that there is no entry for a field
            if !seen_field {
                log::warn!("Entry {} had no values for {field}", self.id())
            }
        }

        for field in E::HierarchicalField::iter() {
            let name = field.to_string();
            let (root, children) = self.tree_field_values(&field);
            entry.add_tree(&name, root, children.iter());
        }

        // All entries should have at least one field
        if !seen_any_field {
            log::error!("Entry {} had no field values at all", self.id());
        }

        for xref in self.cross_references() {
            if xref.name().is_empty() || xref.external_id().is_empty() {
                continue;
            }
            entry.add_ref(xref.name(), xref.external_id())
        }

        // All entries must have a taxid, so we cross reference it correctly
        entry.add_ref("ncbi_taxonomy_id", self.taxid());

        entry
    }
}

#[derive(Debug, Serialize, PartialEq)]
pub struct ValueOnly {
    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct HierarchicalField {
    name: String,
    root: ValueOnly,
    #[serde(rename = "child")]
    children: Vec<ValueOnly>,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct NamedField {
    name: String,

    #[serde(rename = "$value")]
    value: String,
}

#[derive(Debug, Serialize, PartialEq, Default)]
pub struct AdditionalFields {
    #[serde(rename = "field", default)]
    fields: Vec<NamedField>,

    #[serde(rename = "hierarchical_field", default)]
    trees: Vec<HierarchicalField>,
}

#[derive(Debug, Serialize, PartialEq)]
pub struct Ref {
    dbname: String,
    dbkey: String,
}

#[derive(Debug, Serialize, PartialEq, Default)]
pub struct CrossReferences {
    #[serde(rename = "ref", default)]
    refs: Vec<Ref>,
}

#[derive(Debug, Serialize, PartialEq)]
#[serde(rename = "entry")]
pub struct Entry {
    id: String,

    #[serde(rename = "name")]
    name: ValueOnly,

    #[serde(rename = "description")]
    description: ValueOnly,
    cross_references: CrossReferences,
    additional_fields: AdditionalFields,
}

impl<S> From<S> for ValueOnly
where
    S: AsSearchValue,
{
    fn from(value: S) -> Self {
        Self {
            value: value.search_value(),
        }
    }
}

impl Entry {
    pub fn new<N, D>(id: String, name: N, description: D) -> Self
    where
        N: AsSearchValue,
        D: AsSearchValue,
    {
        Self {
            id,
            name: ValueOnly::from(name),
            description: ValueOnly::from(description),
            cross_references: CrossReferences::default(),
            additional_fields: AdditionalFields::default(),
        }
    }

    pub fn add_ref<D, K>(&mut self, database_name: D, key: K)
    where
        D: Into<String>,
        K: AsSearchValue,
    {
        self.cross_references.add_ref(database_name, key);
    }

    pub fn add_field<N, S>(&mut self, name: N, value: S)
    where
        N: Into<String>,
        S: AsSearchValue,
    {
        self.additional_fields.add_field(name, value);
    }

    pub fn add_field_values<N, C, V>(&mut self, name: N, values: V)
    where
        N: Into<String>,
        C: AsSearchValue,
        V: Iterator<Item = C>,
    {
        let name: String = name.into();
        for value in values {
            self.add_field(&name, value);
        }
    }

    pub fn add_tree<N, R, C, I>(&mut self, name: N, root: R, children: I)
    where
        N: Into<String>,
        R: AsSearchValue,
        C: AsSearchValue,
        I: Iterator<Item = C>,
    {
        self.additional_fields.add_tree(name, root, children);
    }
}

impl CrossReferences {
    pub fn add_ref<D, K>(&mut self, database_name: D, key: K)
    where
        D: Into<String>,
        K: AsSearchValue,
    {
        let name: String = database_name.into().replace(" ", "_");
        self.refs.push(Ref {
            dbname: name,
            dbkey: key.search_value(),
        });
    }
}

impl AdditionalFields {
    pub fn add_field<N, V>(&mut self, name: N, value: V)
    where
        N: Into<String>,
        V: AsSearchValue,
    {
        self.fields.push(NamedField {
            name: name.into(),
            value: value.search_value(),
        });
    }

    pub fn add_tree<N, R, C, I>(&mut self, name: N, root: R, children: I)
    where
        N: Into<String>,
        R: AsSearchValue,
        C: AsSearchValue,
        I: Iterator<Item = C>,
    {
        let mut members = Vec::new();
        for child in children {
            members.push(ValueOnly::from(child));
        }
        self.trees.push(HierarchicalField {
            name: name.into(),
            root: root.into(),
            children: members,
        });
    }
}
