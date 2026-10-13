> **Historical record.** These runs predate the current code: the test argument was then
> `log2(domain size)`, and is now the number of validators (use `2^k - 1` for a domain of
> `2^k`), and the test targets are now `basic_377`, `packed_377` and `counting_377`. Output labels have also changed since. Kept for the timings, not as instructions.

Run 12/02/2022 on d26d87767050cb0817481e013cc2f6b8e82264b1

### The basic scheme

#### log(domain_size) = 10, single-threaded
 
> cargo test --release --features "print-trace" --test basic 10

```
Start:   setup
··Start:   Computing 3070 scalars powers
··End:     Computing 3070 scalars powers ...........................................376.300µs
··Start:   3070-scalar mul in G1
··End:     3070-scalar mul in G1 ...................................................409.644ms
··Start:   2-scalar mul in G1
··End:     2-scalar mul in G1 ......................................................8.587ms
End:     setup .....................................................................424.944ms
Start:   signer set commitment
End:     signer set commitment .....................................................109.510ms
Start:   prover precomputation
End:     prover precomputation .....................................................3.092ms
Start:   BW6 prove
End:     BW6 prove .................................................................422.145ms
Start:   BW6 verify
··Start:   linear accountability check
··End:     linear accountability check .............................................99.900µs
··Start:   KZG check
····Start:   linearization polynomial commitment
····End:     linearization polynomial commitment ...................................1.795ms
····Start:   aggregate evaluation claims in zeta
····End:     aggregate evaluation claims in zeta ...................................875.500µs
····Start:   batched KZG openning
····End:     batched KZG openning ..................................................8.737ms
····Start:   lazy subgroup check
····End:     lazy subgroup check ...................................................633.300µs
··End:     KZG check ...............................................................12.882ms
End:     BW6 verify ................................................................14.992ms
```

#### log(domain_size) = 10, multi-threaded

> cargo test --release --features "parallel print-trace" --test basic 10

```
Start:   setup
··Start:   Computing 3070 scalars powers
··End:     Computing 3070 scalars powers ...........................................377.700µs
··Start:   3070-scalar mul in G1
··End:     3070-scalar mul in G1 ...................................................172.988ms
··Start:   2-scalar mul in G1
··End:     2-scalar mul in G1 ......................................................25.782ms
End:     setup .....................................................................206.306ms
Start:   signer set commitment
End:     signer set commitment .....................................................40.431ms
Start:   prover precomputation
End:     prover precomputation .....................................................2.445ms
Start:   BW6 prove
End:     BW6 prove .................................................................159.899ms
Start:   BW6 verify
··Start:   linear accountability check
··End:     linear accountability check .............................................99.400µs
··Start:   KZG check
····Start:   linearization polynomial commitment
····End:     linearization polynomial commitment ...................................2.970ms
····Start:   aggregate evaluation claims in zeta
····End:     aggregate evaluation claims in zeta ...................................1.451ms
····Start:   batched KZG openning
····End:     batched KZG openning ..................................................14.423ms
····Start:   lazy subgroup check
····End:     lazy subgroup check ...................................................1.042ms
··End:     KZG check ...............................................................20.469ms
End:     BW6 verify ................................................................20.995ms
```

#### log(domain_size) = 16, single-threaded

> cargo test --release --features "print-trace" --test basic 16
 
```
Start:   setup
··Start:   Computing 196606 scalars powers
··End:     Computing 196606 scalars powers .........................................13.021ms
··Start:   196606-scalar mul in G1
··End:     196606-scalar mul in G1 .................................................14.938s
··Start:   2-scalar mul in G1
··End:     2-scalar mul in G1 ......................................................9.447ms
End:     setup .....................................................................14.965s
Start:   signer set commitment
End:     signer set commitment .....................................................130.885ms
Start:   prover precomputation
End:     prover precomputation .....................................................21.555ms
Start:   BW6 prove
End:     BW6 prove .................................................................456.569ms
Start:   BW6 verify
··Start:   linear accountability check
··End:     linear accountability check .............................................51.800µs
··Start:   KZG check
····Start:   linearization polynomial commitment
····End:     linearization polynomial commitment ...................................1.848ms
····Start:   aggregate evaluation claims in zeta
····End:     aggregate evaluation claims in zeta ...................................959.700µs
····Start:   batched KZG openning
····End:     batched KZG openning ..................................................9.454ms
····Start:   lazy subgroup check
····End:     lazy subgroup check ...................................................658.700µs
··End:     KZG check ...............................................................14.531ms
End:     BW6 verify ................................................................15.334ms
```

#### log(domain_size) = 16, multi-threaded

> cargo test --release --features "parallel print-trace" --test basic 16

```
Start:   setup
··Start:   Computing 196606 scalars powers
··End:     Computing 196606 scalars powers .........................................21.986ms
··Start:   196606-scalar mul in G1
··End:     196606-scalar mul in G1 .................................................6.166s
··Start:   2-scalar mul in G1
··End:     2-scalar mul in G1 ......................................................16.043ms
End:     setup .....................................................................6.212s
Start:   signer set commitment
End:     signer set commitment .....................................................64.675ms
Start:   prover precomputation
End:     prover precomputation .....................................................28.511ms
Start:   BW6 prove
End:     BW6 prove .................................................................186.508ms
Start:   BW6 verify
··Start:   linear accountability check
··End:     linear accountability check .............................................128.500µs
··Start:   KZG check
····Start:   linearization polynomial commitment
····End:     linearization polynomial commitment ...................................2.862ms
····Start:   aggregate evaluation claims in zeta
····End:     aggregate evaluation claims in zeta ...................................1.377ms
····Start:   batched KZG openning
····End:     batched KZG openning ..................................................13.873ms
····Start:   lazy subgroup check
····End:     lazy subgroup check ...................................................1.008ms
··End:     KZG check ...............................................................19.686ms
End:     BW6 verify ................................................................20.284ms
```