mod helper;


fn main() {
    let validators = helper::parse_args_or(252, "counting_381");
    println!(
        "Running test for the 'counting' scheme on APK-381 (BLS12-381 / BW6-767) for {} validators",
        validators
    );
    apk_proofs::test_helpers::test_counting_scheme_381(validators);
}
