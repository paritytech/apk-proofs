# Fully Succinct BLS Signature Aggregation

Individual BLS signatures on the same message can be aggregated into a single signature
that can be verified in constant time, given the verifier knows the aggregate public key of the set of the actual signers [[1]](https://eprint.iacr.org/2018/483). 
However, computing the aggregate public key is linear in the number of the actual signers and requires the verifier to know the individual public keys.

We avoid such heavy computation for verifiers that are constrained resource-wise and computation-wise (e.g., mobile phones, smart contracts on  blockchains) by designing custom non-interactive succinct arguments of knowledge (SNARKs) that compute and ensure the correctness of an apk, i.e., an aggregated public key of actual signers. This repo contains PoC implementations as well as formalisations for our custom SNARKs for apk, given the verifier knows only a commitment to the list of public keys of all the eligible signers and a bitmask identifying the actual signers of a message.

See the light-client simulation in [bw6/examples](bw6/examples/README.md) for a sketch of a blockchain light client design exploiting such proofs.

Two configurations are implemented: BLS12-377 signatures with proofs over BW6-761 (`Apk377`), and BLS12-381 signatures with proofs over BW6-767 (`Apk381`). The 'packed' scheme is available on the first only.

# Formal Write-up
The formal description and security model for our custom succinct arguments as well as their application to accountable light clients for PoS blockchains can be found [here/stable version](https://eprint.iacr.org/2022/1205) and [here/on-going updates](https://github.com/w3f/apk-proofs/blob/main/Light%20Client.pdf). A high-level summary of this work as well as its connection with related reseach effort from the Web3 Foundation can be found [here](https://research.web3.foundation/en/latest/polkadot/LightClientsBridges/index.html).

# Video Presentations
A video presentation of this work at sub0 2022 is [available here](https://www.youtube.com/watch?v=MCvX9ZZhO4I&list=PLOyWqupZ-WGvywLqJDsMIYdCn8QEa2ShQ&index=19) ([slides](https://docs.google.com/presentation/d/16LlsXWY2Q6_6QGZxkg84evaJqWNk6szX)), and at ZK Summit 7 is [available here](https://www.youtube.com/watch?v=UaPdDYarKGY&list=PLj80z0cJm8QFnY6VLVa84nr-21DNvjWH7&index=19).    

# How to Reproduce Results in the Formal Write-up

1. [Install](https://www.rust-lang.org/tools/install) Rust toolchain.
2. Run one of the commands below from the `bw6` directory.

The argument is the number of validators. The write-up's results are for domains of size 2^10, 2^16 and 2^20, which on BLS12-377/BW6-761 means 2^k - 1 validators (one row of the domain is reserved for the accumulator's initial value).

#### Basic Accountable Scheme 
> cargo test --release --features "parallel print-trace" --test basic_377 1023

> cargo test --release --features "parallel print-trace" --test basic_377 65535

> cargo test --release --features "parallel print-trace" --test basic_377 1048575
 
<br/>

#### Packed Accountable Scheme
> cargo test --release --features "parallel print-trace" --test packed_377 1023

> cargo test --release --features "parallel print-trace" --test packed_377 65535

> cargo test --release --features "parallel print-trace" --test packed_377 1048575

<br/>

#### Counting Scheme
> cargo test --release --features "parallel print-trace" --test counting_377 1023

> cargo test --release --features "parallel print-trace" --test counting_377 65535

> cargo test --release --features "parallel print-trace" --test counting_377 1048575

<br/>

#### The Same Schemes on BLS12-381 / BW6-767
The commands above run on BLS12-377 / BW6-761 (APK-377), which is what the write-up measures. The basic and counting schemes also run on BLS12-381 / BW6-767 (APK-381); the packed scheme does not, as it needs domain sizes divisible by 256 and BW6-767 domain has no such size.

> cargo test --release --features "parallel print-trace" --test basic_381 1023

> cargo test --release --features "parallel print-trace" --test counting_381 1023

APK-381 domains are not powers of two: the validator count plus one rounds up to the next entry of `APK381_DOMAIN_SIZES` (1023 validators use a domain of 1081 = 23 * 47), so the same validator count is not the same domain size on the two curves.

<br/>

The output should look, for example, like this (single-threaded, `--test basic_377 1023`; timings are machine-dependent)
```
Running test for the 'basic' scheme on APK-377 (BLS12-377 / BW6-761) for 1023 validators
Start:   signer set commitment
End:     signer set commitment .....................................................199.144ms
Start:   prover precomputation
End:     prover precomputation .....................................................5.286ms
Start:   prove
End:     prove .....................................................................479.968ms
Start:   verify
··Start:   subgroup checks
··End:     subgroup checks .........................................................835.250µs
··Start:   linear accountability check
··End:     linear accountability check .............................................109.708µs
··Start:   PCS verification
····Start:   linearization polynomial commitment
····End:     linearization polynomial commitment ...................................990.291µs
····Start:   aggregate evaluation claims at zeta
····End:     aggregate evaluation claims at zeta ...................................473.167µs
····Start:   batched PCS opening verification
····End:     batched PCS opening verification ......................................4.876ms
··End:     PCS verification ........................................................6.348ms
End:     verify ....................................................................7.356ms
```

# Running with Docker

The [Dockerfile](Dockerfile) builds the crate and then runs the full test suite and the light-client simulation on both configurations (`apk_377` and `apk_381`, see [bw6/examples](bw6/examples/README.md)). Run from the repository root:

> docker build -t apk-proofs .

> docker run --rm apk-proofs

`docker run --rm apk-proofs tests` runs the tests only, and `docker run --rm apk-proofs light-clients` the light clients only. The light-client size is set with `-e VALIDATORS=1500 -e ERAS=3` (defaults: 1000 validators, 2 eras).
