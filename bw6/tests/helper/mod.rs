use std::env;

/// Parses a validator-set size from the command line.
///
/// Deliberately not a `log_n`: what the caller has is a number of validators, and which
/// evaluation domain that implies is the configuration's business. On APK-381 no domain is a
/// power of two, so a logarithm would name nothing.
///
/// `test_name` is the test target, used only in the hint printed when the argument is missing;
/// which configuration runs is fixed by the calling binary.
pub fn parse_args_or(default_validators: usize, test_name: &str) -> usize {
    match env::args().nth(1) {
        None => {
            println!(
                "VALIDATORS parameter is not provided, using {}. To choose the validator set \
                 size, pass it after the test name, e.g.\n  \
                 cargo test --release --features \"parallel print-trace\" --test {} 1023",
                default_validators, test_name
            );
            default_validators
        }
        Some(arg) => arg
            .parse()
            .unwrap_or_else(|_| panic!("{} is not a valid parameter", arg)),
    }
}
