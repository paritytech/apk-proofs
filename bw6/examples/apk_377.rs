//! The light-client simulation on APK-377: BLS12-377 signatures, proofs over BW6-761, radix-2
//! evaluation domains.
//!
//! Run with `cargo run --release --example apk_377 -- VALIDATORS N_ERAS`. `VALIDATORS` is the
//! size of the validator set, not a logarithm: the configuration works out which evaluation
//! domain that needs. Here it rounds up to a power of two.
//!
//! The body of the simulation is in `examples/common/light_client.rs` and is shared verbatim
//! with `apk_381`. The only difference between the two examples is the config type.

#[path = "common/light_client.rs"]
mod light_client;

use apk_proofs::Bls12_377Config;

fn main() {
    let (validators, eras) = light_client::parse_args("apk_377", 255, 10);
    let rng = &mut ark_std::test_rng(); // Don't use in production code!
    light_client::run::<Bls12_377Config, _>(validators, eras, rng);
}
