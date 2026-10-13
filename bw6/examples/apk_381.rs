//! The light-client simulation on APK-381: BLS12-381 signatures, proofs over BW6-767,
//! mixed-radix evaluation domains.
//!
//! Run with `cargo run --release --example apk_381 -- VALIDATORS N_ERAS`. `VALIDATORS` is the
//! size of the validator set. It is deliberately *not* a `log_n`: BW6-767's scalar field has
//! two-adicity 1, so no domain here is a power of two and a logarithm would name nothing. The
//! configuration picks the domain — a divisor of `3 * 11 * 23 * 47 * 10177`, expanded to `2n`
//! and `6n` — and prints which one it chose.
//!
//! The body of the simulation is in `examples/common/light_client.rs` and is shared verbatim
//! with `apk_377`. The only difference between the two examples is the config type.

#[path = "common/light_client.rs"]
mod light_client;

use apk_proofs::Bls12_381Config;

fn main() {
    let (validators, eras) = light_client::parse_args("apk_381", 252, 10);
    let rng = &mut ark_std::test_rng(); // Don't use in production code!
    light_client::run::<Bls12_381Config, _>(validators, eras, rng);
}
