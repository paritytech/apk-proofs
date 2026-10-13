mod helper;

fn main() {
    let validators = helper::parse_args_or(255, "counting_377");
    println!(
        "Running test for the 'counting' scheme on APK-377 (BLS12-377 / BW6-761) for {} validators",
        validators
    );
    apk_proofs::test_helpers::test_counting_scheme_377(validators);
}
