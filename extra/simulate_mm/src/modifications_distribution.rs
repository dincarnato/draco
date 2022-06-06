use std::fmt;

pub const PROBABILITY: f64 = 0.01927;

pub struct ModificationDistribution {
    n: u32,
    probability: f64,
}

impl ModificationDistribution {
    pub fn new(modificable_bases: u32, probability: f64) -> Self {
        Self {
            n: modificable_bases,
            probability,
        }
    }

    pub fn probability(&self, n_modifications: u32) -> Result<f64, BinomialError> {
        let binom: f64 = binomial_coefficient(self.n, n_modifications)?.into();

        let probability = binom
            * self.probability.powi(n_modifications as i32)
            * (1. - self.probability).powi(self.n.checked_sub(n_modifications).unwrap() as i32);
        Ok(probability)
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
enum BinomialResult {
    Small(u64),
    Big(f64),
}

impl From<BinomialResult> for f64 {
    fn from(binomial_result: BinomialResult) -> Self {
        use BinomialResult::*;
        match binomial_result {
            Small(n) => n as f64,
            Big(n) => n,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BinomialError {
    InvalidSuccesses,
}

impl fmt::Display for BinomialError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("invalid number of successes")
    }
}

fn binomial_coefficient(n: u32, k: u32) -> Result<BinomialResult, BinomialError> {
    use std::collections::BTreeMap;

    if k > n {
        return Err(BinomialError::InvalidSuccesses);
    }

    let mut k_factor_divisors = BTreeMap::new();
    (2..=u64::from(k)).for_each(|mut k| {
        let mut divisor = 2;
        while k > 1 {
            debug_assert!(divisor <= k);
            if k % divisor == 0 {
                *k_factor_divisors.entry(divisor).or_insert(0u64) += 1;
                k /= divisor;
            } else {
                divisor += 1;
            }
        }
    });

    use BinomialResult::*;
    let result = (u64::from(n - k + 1)..=u64::from(n)).fold(Small(1u64), |factor, mut n| {
        for (divisor, occurrences) in &mut k_factor_divisors {
            if n == 1 {
                break;
            }

            while *occurrences > 0 {
                if n % divisor == 0 {
                    n /= divisor;
                    *occurrences -= 1;
                } else {
                    break;
                }
            }
        }
        match factor {
            Small(factor) => factor
                .checked_mul(n)
                .map(Small)
                .unwrap_or_else(|| Big(factor as f64 * n as f64)),
            Big(factor) => Big(factor * n as f64),
        }
    });
    Ok(result)
}
