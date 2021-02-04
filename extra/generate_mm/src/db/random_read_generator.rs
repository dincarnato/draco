use super::Entry;
use bstr::BStr;
use rand::rngs::ThreadRng;

#[derive(Debug)]
pub struct RandomReadGenerator<'a> {
    sequence: &'a BStr,
    modificable_indices: Vec<Vec<u32>>,
    modifications_cdf: Vec<f64>,
    cum_fractions: Vec<f32>,
    rng: ThreadRng,
}

impl<'a> RandomReadGenerator<'a> {
    pub fn new(db_entry: &'a Entry, mut fractions: Vec<f32>, probability: f64) -> Self {
        use crate::modifications_distribution::ModificationDistribution;
        use rand::thread_rng;

        fractions
            .iter_mut()
            .scan(0f32, |cum, cur| {
                *cum += *cur;
                Some((cur, *cum))
            })
            .for_each(|(cur, cum)| *cur = cum);
        assert!((fractions.last().expect("fractions vector cannot be empty") - 1.).abs() <= 0.001);
        fractions.truncate(fractions.len() - 1);

        let mod_dist = ModificationDistribution::new(db_entry.sequence.len() as u32, probability);
        let modifications_cdf = (0u32..)
            .map(|k| mod_dist.probability(k).unwrap())
            .scan(0., |cdf, p| {
                *cdf += p;
                Some(*cdf)
            })
            .take_while(|&cdf| cdf <= 0.9999)
            .collect();

        let modificable_indices: Vec<_> = db_entry
            .profiles
            .iter()
            .map(|profile| &profile.0)
            .map(|profile| {
                profile
                    .iter()
                    .enumerate()
                    .filter(|(_, &modificability)| modificability > 0)
                    .map(|(index, _)| index as u32)
                    .collect::<Vec<_>>()
            })
            .collect();
        assert!(modificable_indices
            .iter()
            .all(
                |modificable_indices| modificable_indices.iter().all(|&index| {
                    let base = db_entry.sequence[index as usize];
                    base == b'A' || base == b'C'
                })
            ));
        assert_ne!(modificable_indices.len(), 0);
        assert_eq!(modificable_indices.len() - 1, fractions.len());

        Self {
            sequence: &db_entry.sequence,
            modificable_indices,
            modifications_cdf,
            cum_fractions: fractions,
            rng: thread_rng(),
        }
    }
}

impl Iterator for RandomReadGenerator<'_> {
    type Item = (usize, Vec<u32>);

    fn next(&mut self) -> Option<Self::Item> {
        use rand::{seq::SliceRandom, Rng};

        let random_fraction = self.rng.gen();
        let profile_index = self
            .cum_fractions
            .iter()
            .position(|&fraction| fraction >= random_fraction)
            .unwrap_or_else(|| self.modificable_indices.len() - 1);
        let modificable_indices = &self.modificable_indices[profile_index];

        let p: f64 = self.rng.gen();
        let n_modifications = self
            .modifications_cdf
            .iter()
            .position(|&cdf| cdf > p)
            .unwrap_or_else(|| self.modifications_cdf.len());

        let mut sampled: Vec<_> = modificable_indices
            .choose_multiple(&mut self.rng, n_modifications)
            .cloned()
            .collect();
        sampled.sort_unstable();
        Some((profile_index, sampled))
    }
}
