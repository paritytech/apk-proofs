## Light-client simulation

A chain whose validator set rotates every era, and a light client that follows it by checking
one APK proof and one aggregate BLS signature per era. The design is explained in the comments
of [common/light_client.rs](common/light_client.rs), which is written once against `ApkConfig`
and run by two examples that differ only in the configuration:

- [apk_377.rs](apk_377.rs): BLS12-377 signatures, proofs over BW6-761, radix-2 domains.
- [apk_381.rs](apk_381.rs): BLS12-381 signatures, proofs over BW6-767, mixed-radix domains.

Run from the `bw6` directory. `VALIDATORS` is the size of the validator set;
the configuration chooses the evaluation domain and prints it.

> cargo run --release --features "parallel print-trace" --example apk_381 -- VALIDATORS N_ERAS

For example, `--example apk_381 -- 252 1` (single-threaded; timings are machine-dependent)
produces

```
The configuration picked an evaluation domain of 253 points for 252 validators.

Start:   Generating PCS params to support 252 signers
End:     Generating PCS params to support 252 signers ..............................96.886ms

Genesis: validator set size = 252, quorum = 169

Start:   Computing commitment to the set of initial 252 validators
End:     Computing commitment to the set of initial 252 validators .................69.400ms

Era 1

Start:   Each (honest) validator computes the commitment to the new validator set of size 252 and signs the commitment
End:     Each (honest) validator computes the commitment to the new validator set of size 252 and signs the commitment 255.544ms

Start:   Helper aggregates 200 individual signatures on the same commitment and generates accountable light client proof
End:     Helper aggregates 200 individual signatures on the same commitment and generates accountable light client proof 285.193ms

Start:   Light client verifies light client proof for 200 signers
··Start:   apk proof verification
····Start:   subgroup checks
····End:     subgroup checks .......................................................2.653ms
····Start:   linear accountability check
····End:     linear accountability check ...........................................42.167µs
····Start:   PCS verification
······Start:   linearization polynomial commitment
······End:     linearization polynomial commitment .................................1.035ms
······Start:   aggregate evaluation claims at zeta
······End:     aggregate evaluation claims at zeta .................................505.916µs
······Start:   batched PCS opening verification
······End:     batched PCS opening verification ....................................5.301ms
····End:     PCS verification ......................................................6.851ms
··End:     apk proof verification ..................................................10.782ms
··Start:   aggregate BLS signature verification
··End:     aggregate BLS signature verification ....................................1.585ms
End:     Light client verifies light client proof for 200 signers ..................12.374ms

Era 2

Start:   Each (honest) validator computes the commitment to the new validator set of size 252 and signs the commitment
End:     Each (honest) validator computes the commitment to the new validator set of size 252 and signs the commitment 259.974ms

Start:   Helper aggregates 207 individual signatures on the same commitment and generates accountable light client proof
End:     Helper aggregates 207 individual signatures on the same commitment and generates accountable light client proof 278.119ms

Start:   Light client verifies light client proof for 207 signers
··Start:   apk proof verification
····Start:   subgroup checks
····End:     subgroup checks .......................................................2.659ms
····Start:   linear accountability check
····End:     linear accountability check ...........................................38.708µs
····Start:   PCS verification
······Start:   linearization polynomial commitment
······End:     linearization polynomial commitment .................................1.012ms
······Start:   aggregate evaluation claims at zeta
······End:     aggregate evaluation claims at zeta .................................484.083µs
······Start:   batched PCS opening verification
······End:     batched PCS opening verification ....................................5.256ms
····End:     PCS verification ......................................................6.766ms
··End:     apk proof verification ..................................................10.678ms
··Start:   aggregate BLS signature verification
··End:     aggregate BLS signature verification ....................................1.516ms
End:     Light client verifies light client proof for 207 signers ..................12.202ms
```

The parameters are generated locally from a known trapdoor (`InsecureSetup`), which is fine for a
simulation and unsafe anywhere else.
