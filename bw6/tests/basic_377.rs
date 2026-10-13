mod helper;

fn main() {
    let validators = helper::parse_args_or(255, "basic_377");
    println!(
        "Running test for the 'basic' scheme on APK-377 (BLS12-377 / BW6-761) for {} validators",
        validators
    );
    apk_proofs::test_helpers::test_simple_scheme_377(validators);
}
