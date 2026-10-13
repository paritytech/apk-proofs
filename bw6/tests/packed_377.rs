mod helper;

fn main() {
    let validators = helper::parse_args_or(255, "packed_377");
    println!(
        "Running test for the 'packed' scheme on APK-377 (BLS12-377 / BW6-761) for {} validators",
        validators
    );
    apk_proofs::test_helpers::test_packed_scheme_377(validators);
}
