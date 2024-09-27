mod db;
mod modifications_distribution;
mod mutation_map;

use std::{
    error::Error,
    fmt::{self, Display},
    num::ParseFloatError,
    path::PathBuf,
    str::FromStr,
};
use structopt::{clap::ArgGroup, StructOpt};

#[derive(StructOpt)]
#[structopt(group = ArgGroup::with_name("reads_amount").required(true))]
struct ArgsOpt {
    /// Input file with structure profiles
    #[structopt(parse(from_os_str))]
    db_file: PathBuf,

    /// Output file with structure profiles updated according to the simulation
    #[structopt(parse(from_os_str), short = "o", long = "outProfiles")]
    db_out: PathBuf,

    /// Output mutation map (MM) file
    #[structopt(parse(from_os_str))]
    mm_file: PathBuf,

    /// Number of reads mapping to each transcript
    /// [Note: this parameter and "--meanCoverage" are mutually exclusive]
    #[structopt(short = "n", long = "nReads", group = "reads_amount")]
    n_reads: Option<usize>,

    /// Mean sequencing depth (coverage) per base
    /// [Note: this parameter and "--nReads" are mutually exclusive]
    #[structopt(short = "c", long = "meanCoverage", group = "reads_amount")]
    mean_coverage: Option<u32>,

    /// Length (in bp) of the simulated reads
    #[structopt(short = "s", long = "readLen")]
    read_size: u32,

    /// Comma-separated list of % conformation stoichiometries
    /// [Note: the stoichiometries must sum to approx. 100 (tollerance: 97-103).
    /// When no stoichiometry is specified, the conformations are assumed to be equimolar]
    #[structopt(
        short = "p",
        long = "stoichiometry",
        parse(try_from_str = parse_percentages)
    )]
    fractions: Vec<Vec<f32>>,

    /// Output MM file "human-readable" version
    #[structopt(short, long)]
    text: Option<PathBuf>,

    /// Sets the `p` value for generation of the binomial distribution of mutations
    /// [Note: the default value (0.01927) has been learnt empirically from Homan et al., 2014]
    #[structopt(default_value, long)]
    probability: Probability,
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
struct Probability(f64);

impl Default for Probability {
    fn default() -> Self {
        Self(modifications_distribution::PROBABILITY)
    }
}

impl fmt::Display for Probability {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl FromStr for Probability {
    type Err = <f64 as FromStr>::Err;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        f64::from_str(s).map(Self)
    }
}

#[derive(Debug, Clone)]
enum ParsePercentageError {
    InvalidFloat(ParseFloatError),
    InvalidSum,
}

impl Display for ParsePercentageError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let s = match self {
            ParsePercentageError::InvalidFloat(_) => "percentage is an invalid float",
            ParsePercentageError::InvalidSum => "percentages do not sum to 100",
        };
        f.write_str(s)
    }
}

impl Error for ParsePercentageError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ParsePercentageError::InvalidFloat(source) => Some(source),
            ParsePercentageError::InvalidSum => None,
        }
    }
}

fn parse_percentages(src: &str) -> Result<Vec<f32>, ParsePercentageError> {
    let mut fractions = src
        .split(',')
        .map(|raw| {
            raw.parse()
                .map(|percentage: f32| percentage / 100.)
                .map_err(ParsePercentageError::InvalidFloat)
        })
        .collect::<Result<Vec<_>, _>>()?;
    let sum: f32 = fractions.iter().sum();

    if (sum - 1.).abs() > 0.03 {
        Err(ParsePercentageError::InvalidSum)
    } else {
        fractions.iter_mut().for_each(|fraction| *fraction /= sum);
        Ok(fractions)
    }
}

fn main() {
    use db::Db;
    use rand::{thread_rng, Rng};
    use std::{
        fs::File,
        io::{BufWriter, Write},
        ops::Range,
    };

    let args_opt = ArgsOpt::from_args();
    let db = Db::new(&args_opt.db_file).expect("cannot read DB file");
    let mut mm_file =
        BufWriter::new(File::create(&args_opt.mm_file).expect("cannot create MM file"));
    let mut rng = thread_rng();

    let mut db_out =
        BufWriter::new(File::create(&args_opt.db_out).expect("cannot create the output DB"));

    let mut text_file = args_opt
        .text
        .as_ref()
        .map(|text| File::create(text).expect("cannot create the output text file"));

    let probability = args_opt.probability.0;

    db.entries
        .into_iter()
        .enumerate()
        .for_each(|(entry_index, entry)| {
            let range_max = entry.sequence.len() as u32 - args_opt.read_size;
            let mut modifications = vec![vec![0; entry.sequence.len()]; entry.profiles.len()];

            let len_profiles = entry.profiles.len();
            let fractions = args_opt
                .fractions
                .iter()
                .find(|fractions| fractions.len() == len_profiles)
                .cloned()
                .unwrap_or_else(|| vec![1. / len_profiles as f32; len_profiles]);

            let profile_mapper = |(profile_index, mut modification_indices, range): (
                usize,
                Vec<u32>,
                Range<u32>,
            )| {
                modification_indices.retain(|&index| index >= range.start && index < range.end);
                let cur_modifications = &mut modifications[profile_index];
                modification_indices
                    .iter()
                    .for_each(|&base_index| cur_modifications[base_index as usize] += 1);

                mutation_map::Read::new(range.start, range.end, modification_indices)
            };

            let mutations_generator = entry
                .random_read_generator(fractions, probability)
                .map(|(profile_index, mutations_indices)| {
                    let begin = if range_max == 0 {
                        0
                    } else {
                        rng.gen_range(0..range_max)
                    };
                    let end = begin + args_opt.read_size;
                    (profile_index, mutations_indices, begin..end)
                });

            let reads: Vec<_> = {
                if let Some(n_reads) = args_opt.n_reads {
                    mutations_generator.take(n_reads).map(profile_mapper).collect()
                } else {
                    let mean_coverage = args_opt.mean_coverage.unwrap();

                    let initial_data = (false, vec![0u32; entry.sequence.len()]);
                    mutations_generator.scan(initial_data, |(minimally_covered, coverages), (profile_index, modification_indices, range)| {
                        if *minimally_covered {
                            let mean = coverages.iter().map(|&coverage| coverage as f32 / entry.sequence.len() as f32).sum::<f32>() as u32;
                            if mean >= mean_coverage {
                                return None
                            }
                        }

                        let coverages_range = &mut coverages[range.start as usize..range.end as usize];
                        coverages_range.iter_mut().for_each(|coverage| *coverage += 1);
                        if !*minimally_covered && coverages_range.iter().any(|&coverage| coverage >= mean_coverage) {
                            *minimally_covered = true;
                        }

                        Some((profile_index, modification_indices, range))
                    }).map(profile_mapper).collect()
                }
            };

            db_out
                .write_all(entry.sequence.as_slice())
                .expect("cannot write to output DB");
            db_out
                .write_all(b"\n")
                .expect("cannot write to output DB");
            entry.profiles.iter().for_each(|profile| {
                let modificability: Vec<_> = profile
                    .0
                    .iter()
                    .map(|&count| if count > 0 { b'.' } else { b'x' })
                    .chain(std::iter::once(b'\n'))
                    .collect();
                db_out
                    .write_all(&modificability)
                    .expect("cannot write to output DB");
            });
            modifications.into_iter().for_each(|modifications| {
                writeln!(db_out, "{}", itertools::join(modifications, ","))
                    .expect("cannot write to output DB");
            });
            db_out
                .write_all(b"\n")
                .expect("cannot write to output DB");

            let transcript_name =
                format!("RNA_{}", entry_index + 1).into();

            if let Some(text_file) = &mut text_file {
                writeln!(text_file, "{}\n{}", transcript_name, entry.sequence).unwrap();
                for read in &reads {
                    write!(text_file, "{}-{}", read.begin, read.end).unwrap();

                    let mut indices = read.indices.iter();
                    if let Some(index) = indices.next() {
                        write!(text_file, " {}", index).unwrap();
                        for index in indices {
                            write!(text_file, ",{}", index).unwrap();
                        }
                    }
                    writeln!(text_file).unwrap();
                }
            }

            let transcript = mutation_map::Transcript::new(
                transcript_name,
                entry.sequence,
                reads,
            );

            transcript
                .serialize_into(&mut mm_file)
                .expect("cannot write transcript to MM file");
        });

    mutation_map::MutationMap::write_end_marker(mm_file)
        .expect("cannot write end marker of MM file");
}
