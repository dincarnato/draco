mod random_read_generator;

use bstr::BString;
use random_read_generator::RandomReadGenerator;
use std::{io, path::Path, str};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Db {
    pub entries: Vec<Entry>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Entry {
    pub sequence: BString,
    pub profiles: Vec<Profile>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Profile(pub Vec<u32>);

#[derive(Debug, Clone, PartialEq, Eq)]
struct ModificabilityProfile(Vec<bool>);

impl Db {
    pub fn new(filename: &Path) -> io::Result<Self> {
        use bstr::io::BufReadExt;
        use std::{fs::File, io::BufReader};

        let buffer = BufReader::new(File::open(filename)?);
        let mut lines = buffer.byte_lines();
        let mut entries = Vec::new();
        while let Some(sequence) = lines.next() {
            let sequence = sequence?.into();
            let mut modificability_profiles: Vec<ModificabilityProfile> = Vec::new();

            let mut line = loop {
                let line = lines.next().unwrap()?;
                match line[0] {
                    b'(' | b')' | b'x' | b'.' => modificability_profiles.push(line.into()),
                    _ => break line,
                }
            };

            let mut profiles = Vec::new();
            for modificability_profile in modificability_profiles {
                let profile: Profile = line.into();
                assert!(profile
                    .0
                    .iter()
                    .zip(modificability_profile.0)
                    .all(|(&count, modificable)| (count != 0) == modificable));

                profiles.push(profile);
                line = lines.next().unwrap()?;
            }
            assert!(line.is_empty());

            entries.push(Entry { sequence, profiles });
        }

        Ok(Self { entries })
    }
}

impl From<Vec<u8>> for ModificabilityProfile {
    fn from(s: Vec<u8>) -> Self {
        let profile: Vec<_> = s
            .iter()
            .map(|c| match c {
                b'x' => false,
                b'(' => false,
                b')' => false,
                b'.' => true,
                _ => panic!("invalid modificability profile character '{}'", c),
            })
            .collect();

        Self(profile)
    }
}

impl<BS> From<BS> for Profile
where
    BS: AsRef<[u8]>,
{
    fn from(s: BS) -> Self {
        let profile: Vec<u32> = s
            .as_ref()
            .split(|&c| c == b',')
            .map(|c| {
                str::from_utf8(c)
                    .expect("invalid UTF-8 string")
                    .parse()
                    .expect("expected positive integer")
            })
            .collect();
        Self(profile)
    }
}

impl Entry {
    pub fn random_read_generator(
        &self,
        fractions: Vec<f32>,
        probability: f64,
    ) -> RandomReadGenerator {
        RandomReadGenerator::new(self, fractions, probability)
    }
}
