use std::{
    fs::File,
    io,
    io::BufRead,
    path::Path,
};

use quick_xml::{
    events::Event,
    Reader,
};

use fallible_iterator::FallibleIterator;
use thiserror::Error;

use crate::publications::reference::{
    Author,
    AuthorBuilder,
    AuthorBuildingError,
    Reference,
    ReferenceBuildError,
    ReferenceBuilder,
};

#[derive(Error, Debug)]
pub enum XmlError {
    #[error("Cannot handle pmc property {0}")]
    UnknownPmcProperty(String),

    #[error("Cannot handle article property {0}")]
    UnknownArticleProperty(String),

    #[error("Cannot handle article property {0}")]
    UnknownAuthorEvent(String),

    #[error("Could not build authors")]
    AuthorError(#[from] AuthorBuildingError),

    #[error("Could not build the Reference")]
    ReferenceError(#[from] ReferenceBuildError),

    #[error("Could not parse XML")]
    XmlIssue(#[from] quick_xml::Error),

    #[error("IO Error")]
    IoError(#[from] io::Error),
}

pub struct XmlIterator<R: BufRead> {
    reader: Reader<R>,
}

fn node_name(raw: &[u8]) -> String {
    String::from_utf8(raw.to_vec()).unwrap()
}

impl<'a> XmlIterator<&'a [u8]> {
    pub fn from_str(s: &'a str) -> XmlIterator<&'a [u8]> {
        Self {
            reader: Reader::from_reader(s.as_bytes()),
        }
    }
}

impl XmlIterator<io::BufReader<File>> {
    pub fn from_file<P: AsRef<Path>>(
        path: P,
    ) -> Result<XmlIterator<io::BufReader<File>>, XmlError> {
        let file = File::open(path)?;
        let reader = io::BufReader::new(file);
        Ok(Self::from_reader(reader))
    }
}

impl<R: BufRead> XmlIterator<R> {
    pub fn from_reader(reader: R) -> Self {
        Self {
            reader: Reader::from_reader(reader),
        }
    }

    fn next_author(&mut self) -> Result<Option<Author>, XmlError> {
        let mut buf = Vec::new();
        let mut builder = AuthorBuilder::new();
        loop {
            match self.reader.read_event(&mut buf)? {
                Event::End(ref e) => match e.name() {
                    b"Author" => break,
                    b"AuthorList" => {
                        return Ok(None);
                    },
                    b"LastName" => (),
                    b"Initials" => (),
                    n => {
                        return Err(XmlError::UnknownAuthorEvent(node_name(n)));
                    },
                },
                Event::Start(ref e) => match e.name() {
                    b"Author" => (),
                    b"LastName" => {
                        builder.set_last_name(self.reader.read_text(e.name(), &mut Vec::new())?)
                    },
                    b"Initials" => {
                        builder.set_first_name(self.reader.read_text(e.name(), &mut Vec::new())?)
                    },
                    n => {
                        return Err(XmlError::UnknownAuthorEvent(node_name(n)));
                    },
                },
                _ => (),
            }
            buf.clear();
        }

        Ok(Some(builder.build()?))
    }
}

impl<R: BufRead> FallibleIterator for XmlIterator<R> {
    type Error = XmlError;
    type Item = Reference;

    fn next(&mut self) -> Result<Option<Self::Item>, Self::Error> {
        let mut buf = Vec::new();
        loop {
            match self.reader.read_event(&mut buf)? {
                Event::Start(ref e) => match e.name() {
                    b"PMC_ARTICLE" => {
                        let mut builder = ReferenceBuilder::new();
                        loop {
                            match self.reader.read_event(&mut buf)? {
                                Event::Start(ref e) => match e.name() {
                                    b"id" => (),
                                    b"source" => (),
                                    b"issue" => (),
                                    b"journalVolume" => (),
                                    b"journalIssn" => (),
                                    b"pageInfo" => (),
                                    b"pubType" => (),
                                    b"pmid" => builder.set_pmid(
                                        self.reader.read_text(e.name(), &mut Vec::new())?,
                                    ),
                                    b"pmcid" => builder.set_pmcid(
                                        self.reader.read_text(e.name(), &mut Vec::new())?,
                                    ),
                                    b"DOI" => builder
                                        .set_doi(self.reader.read_text(e.name(), &mut Vec::new())?),
                                    b"title" => builder.set_title(
                                        self.reader.read_text(e.name(), &mut Vec::new())?,
                                    ),
                                    b"AuthorList" => {
                                        while let Some(author) = self.next_author()? {
                                            builder.add_author(author);
                                        }
                                    },
                                    b"journalTitle" => builder.set_journal(
                                        self.reader.read_text(e.name(), &mut Vec::new())?,
                                    ),
                                    b"pubYear" => builder.set_year(
                                        self.reader.read_text(e.name(), &mut Vec::new())?,
                                    ),
                                    n => {
                                        return Err(XmlError::UnknownArticleProperty(node_name(n)));
                                    },
                                },
                                Event::End(ref e) => match e.name() {
                                    b"PMC_ARTICLE" => break,
                                    _ => (),
                                },
                                _ => (),
                            }
                        }
                        return Ok(Some(builder.build()?));
                    },
                    b"PMCSet" => (),
                    n => {
                        return Err(XmlError::UnknownPmcProperty(node_name(n)));
                    },
                },
                Event::Eof => return Ok(None),
                _ => (),
            }
            buf.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::publications::reference::Author;

    #[test]
    fn can_parse_simple_xml() -> Result<(), XmlError> {
        let xml = r#"
        <PMCSet>
            <PMC_ARTICLE>
                    <id>26184978</id>
                    <source>MED</source>
                    <pmid>26184978</pmid>
                    <pmcid>PMC4505325</pmcid>
                    <DOI>10.1038/srep12276</DOI>
                    <title>MiR-135b-5p and MiR-499a-3p Promote Cell Proliferation and Migration in Atherosclerosis by Directly Targeting MEF2C.</title>
                    <AuthorList CompleteYN="Y">
                            <Author ValidYN="Y">
                                    <LastName>Xu</LastName>
                                    <Initials>Z</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Han</LastName>
                                    <Initials>Y</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Liu</LastName>
                                    <Initials>J</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Jiang</LastName>
                                    <Initials>F</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Hu</LastName>
                                    <Initials>H</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Wang</LastName>
                                    <Initials>Y</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Liu</LastName>
                                    <Initials>Q</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Gong</LastName>
                                    <Initials>Y</Initials>
                            </Author>
                            <Author ValidYN="Y">
                                    <LastName>Li</LastName>
                                    <Initials>X</Initials>
                            </Author>
                    </AuthorList>
                    <journalTitle>Scientific reports</journalTitle>
                    <issue></issue>
                    <journalVolume>5</journalVolume>
                    <pubYear>2015</pubYear>
                    <journalIssn></journalIssn>
                    <pageInfo>12276</pageInfo>
                    <pubType>&quot;Journal Article&quot;, &quot;Research Support, Non-U.S. Gov&apos;t&quot;</pubType>
            </PMC_ARTICLE>
        </PMCSet>
        "#;

        let mut iterator = XmlIterator::from_str(xml);
        let mut refs = Vec::new();
        while let Some(reference) = iterator.next()? {
            refs.push(reference);
        }

        let mut builder = Reference::builder();
        builder.set_pmid(String::from("26184978"));
        builder.set_pmcid(String::from("PMC4505325"));
        builder.set_doi(String::from("10.1038/srep12276"));
        builder.set_year(String::from("2015"));
        builder.set_journal(String::from("Scientific reports"));
        builder.set_title(String::from("MiR-135b-5p and MiR-499a-3p Promote Cell Proliferation and Migration in Atherosclerosis by Directly Targeting MEF2C."));
        builder.add_author(Author::from(("Z", "Xu")));
        builder.add_author(Author::from(("Y", "Han")));
        builder.add_author(Author::from(("J", "Liu")));
        builder.add_author(Author::from(("F", "Jiang")));
        builder.add_author(Author::from(("H", "Hu")));
        builder.add_author(Author::from(("Y", "Wang")));
        builder.add_author(Author::from(("Q", "Liu")));
        builder.add_author(Author::from(("Y", "Gong")));
        builder.add_author(Author::from(("X", "Li")));

        assert_eq!(refs, vec![builder.build()?]);
        Ok(())
    }
}
