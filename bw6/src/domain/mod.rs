//! Evaluation domains and the FFT strategies behind them.
//!
//! Three layers, from the inside out:
//!
//! - [`FftDomain`] is a single multiplicative subgroup with transforms over it. [`Radix2Domain`]
//!   delegates to arkworks; [`CooleyTukeyDomain`] handles the mixed-radix sizes BW6-767 forces;
//!   `NaiveDomain`, the O(n^2) reference the other two are tested against, exists only under
//!   `test` and the `test-utils` feature.
//! - [`DomainSet`] is the *triple* the protocol actually works over, chosen together because the
//!   three sizes constrain each other. [`Radix2DomainSet`] is `n, 2n, 4n`; [`SmoothDomainSet`] is
//!   `n, 2n, 6n` over a [`DomainSizes`] list precomputed for the field.
//! - Everything above the domain layer takes one `D: DomainSet` parameter and never asks which
//!   curve it is on. See [`crate::config`].
//!
//! See [`FftDomain`] for why the two configurations cannot share one domain implementation.

mod cooley_tukey;
#[cfg(any(test, feature = "test-utils"))]
mod naive;
mod radix2;
mod set;
mod subgroup;
mod types;

pub use cooley_tukey::CooleyTukeyDomain;
#[cfg(any(test, feature = "test-utils"))]
pub use naive::NaiveDomain;
pub use radix2::Radix2Domain;
pub(crate) use set::nesting_index;
pub use set::{DomainSet, DomainSizes, DomainTriple, Radix2DomainSet, SmoothDomainSet};
pub use subgroup::subgroup_generator;
pub use types::{DomainError, FftDomain, SupportsPackedScheme};
