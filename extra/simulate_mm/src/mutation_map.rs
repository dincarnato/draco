use bincode::serialize_into;
use bstr::BString;
use std::{io::Write, os::raw::c_char};

#[derive(Debug)]
pub struct MutationMap {
    transcripts: Vec<Transcript>,
}

#[derive(Debug)]
pub struct Transcript {
    pub id: BString,
    pub sequence: BString,
    pub reads: Vec<Read>,
}

#[derive(Debug)]
pub struct Read {
    pub begin: u32,
    pub end: u32,
    pub indices: Vec<u32>,
}

impl MutationMap {
    #[allow(dead_code)]
    pub fn serialize_into<W>(&self, mut writer: W) -> bincode::Result<()>
    where
        W: Write,
    {
        for transcript in &self.transcripts {
            transcript.serialize_into(&mut writer)?;
        }

        MutationMap::write_end_marker(writer)
    }

    pub fn write_end_marker<W>(mut writer: W) -> bincode::Result<()>
    where
        W: Write,
    {
        serialize_into(&mut writer, &b'\x5b')?;
        serialize_into(&mut writer, &b'\x6d')?;
        serialize_into(&mut writer, &b'\x6d')?;
        serialize_into(&mut writer, &b'\x65')?;
        serialize_into(&mut writer, &b'\x6f')?;
        serialize_into(&mut writer, &b'\x66')?;
        serialize_into(&mut writer, &b'\x5d')
    }
}

impl Transcript {
    pub fn new(id: BString, sequence: BString, reads: Vec<Read>) -> Self {
        Self {
            id,
            sequence,
            reads,
        }
    }

    pub fn serialize_into<W>(&self, mut writer: W) -> bincode::Result<()>
    where
        W: Write,
    {
        serialize_into(&mut writer, &(self.id.len() as u16 + 1))?;
        for &c in self.id.as_slice() {
            serialize_into(&mut writer, &(c as c_char))?;
        }
        serialize_into(&mut writer, &(0 as c_char))?;

        serialize_into(&mut writer, &(self.sequence.len() as u32))?;
        for chunk in self.sequence.as_slice().chunks(2) {
            let c = if chunk.len() == 2 {
                (get_base_repr(chunk[0]) << 4) | get_base_repr(chunk[1])
            } else {
                get_base_repr(chunk[0]) << 4
            };
            serialize_into(&mut writer, &c)?;
        }

        serialize_into(&mut writer, &(self.reads.len() as u32))?;
        for read in &self.reads {
            read.serialize_into(&mut writer)?;
        }
        Ok(())
    }
}

impl Read {
    pub fn new(begin: u32, end: u32, indices: Vec<u32>) -> Self {
        assert!(begin <= end);
        // let size = end - begin;
        // assert!(indices.iter().all(|&index| index < size));

        Self {
            begin,
            end,
            indices,
        }
    }

    fn serialize_into<W>(&self, mut writer: W) -> bincode::Result<()>
    where
        W: Write,
    {
        serialize_into(&mut writer, &self.begin)?;
        serialize_into(&mut writer, &(self.end - 1))?;
        serialize_into(&mut writer, &(self.indices.len() as u32))?;
        for index in &self.indices {
            serialize_into(&mut writer, index as &u32)?;
        }
        Ok(())
    }
}

fn get_base_repr(base: u8) -> u8 {
    match base {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        b'N' => 4,
        _ => panic!("invalid base"),
    }
}
