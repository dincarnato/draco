use std::convert::TryFrom;

use super::Entry;
use rand::{rngs::ThreadRng, seq::SliceRandom, Rng};

#[derive(Debug)]
pub struct RandomReadGenerator<'a> {
    modificable_indices: Vec<Vec<u32>>,
    modifications_cdf: Vec<f64>,
    cum_fractions: Vec<f32>,
    sample_fn: fn(&mut Self, usize, usize) -> Vec<u32>,
    db_entry: &'a Entry,
    rng: ThreadRng,
}

impl<'a> RandomReadGenerator<'a> {
    pub fn new(
        db_entry: &'a Entry,
        mut fractions: Vec<f32>,
        probability: f64,
        profile_weights: bool,
    ) -> Self {
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
        assert_ne!(modificable_indices.len(), 0);
        assert_eq!(modificable_indices.len() - 1, fractions.len());

        let sample_fn = if profile_weights {
            Self::sample_indices_weighted
        } else {
            Self::sample_indices_equally
        };

        Self {
            modificable_indices,
            modifications_cdf,
            cum_fractions: fractions,
            sample_fn,
            db_entry,
            rng: thread_rng(),
        }
    }

    fn sample_indices_equally(&mut self, profile_index: usize, n_modifications: usize) -> Vec<u32> {
        self.modificable_indices[profile_index]
            .choose_multiple(&mut self.rng, n_modifications)
            .cloned()
            .collect()
    }

    fn sample_indices_weighted(
        &mut self,
        profile_index: usize,
        n_modifications: usize,
    ) -> Vec<u32> {
        let profile = &self.db_entry.profiles[profile_index];
        self.modificable_indices[profile_index]
            .choose_multiple_weighted(&mut self.rng, n_modifications, |&index| {
                profile.0[usize::try_from(index).unwrap()]
            })
            .expect("unable to choose multiple indices using profile as weights")
            .cloned()
            .collect()
    }
}

impl Iterator for RandomReadGenerator<'_> {
    type Item = (usize, Vec<u32>);

    fn next(&mut self) -> Option<Self::Item> {
        let random_fraction = self.rng.gen();
        let profile_index = self
            .cum_fractions
            .iter()
            .position(|&fraction| fraction >= random_fraction)
            .unwrap_or(self.modificable_indices.len() - 1);

        let p: f64 = self.rng.gen();
        let n_modifications = self
            .modifications_cdf
            .iter()
            .position(|&cdf| cdf > p)
            .unwrap_or(self.modifications_cdf.len());

        let mut sampled = (self.sample_fn)(self, profile_index, n_modifications);
        sampled.sort_unstable();
        Some((profile_index, sampled))
    }
}
