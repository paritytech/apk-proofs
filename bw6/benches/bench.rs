extern crate apk_proofs;

use apk_proofs::instances::bls12_377_bw6_761::{OMEGA, U};
use apk_proofs::{setup, FftDomain, Keyset};
use ark_bw6_761::{Fr, BW6_761};
use ark_ec::VariableBaseMSM;
use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{FftField, Field};
use ark_poly::univariate::DensePolynomial;
use ark_poly::{DenseUVPolynomial, EvaluationDomain, Radix2EvaluationDomain};
use ark_std::{test_rng, UniformRand};
use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use w3f_pcs::pcs::kzg::KZG;
use w3f_pcs::pcs::PcsParams;

/// BW6-761's scalar field is 2-adic, so the benches use the radix-2 domain.
type BenchDomain = apk_proofs::Radix2DomainSet<Fr>;

fn barycentric_evaluation<F: Field>(c: &mut Criterion, n: u32) {
    use ark_poly::{EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain};

    let rng = &mut test_rng();
    let n = std::convert::TryInto::try_into(n).unwrap();
    let domain = Radix2EvaluationDomain::new(n).unwrap();
    let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();
    let evals2 = evals.clone();
    let z = Fr::rand(rng);

    c.bench_function("barycentric_evaluation", move |b| {
        b.iter(|| {
            apk_proofs::utils::barycentric_eval_at(
                black_box(z),
                black_box(&evals),
                black_box(&apk_proofs::Radix2Domain(domain)),
            )
        })
    });

    let evals = Evaluations::from_vec_and_domain(evals2, domain);
    c.bench_function("interpolate + evaluate", move |b| {
        b.iter(|| {
            black_box(&evals)
                .interpolate_by_ref()
                .evaluate(black_box(&z));
        })
    });
}

fn msm<G: AffineRepr>(c: &mut Criterion, n: usize) {
    let rng = &mut test_rng();

    let nu = G::ScalarField::rand(rng);
    let scalars = (0..n).map(|i| nu.pow([i as u64])).collect::<Vec<_>>();
    let bases = (0..n)
        .map(|_| G::Group::rand(rng).into_affine())
        .collect::<Vec<_>>();

    {
        let (scalars, bases) = (scalars.clone(), bases.clone());
        c.bench_function("ark_ec::msm::VariableBaseMSM", move |b| {
            b.iter(|| <G::Group as VariableBaseMSM>::msm(black_box(&bases), black_box(&scalars)))
        });
    }

    {
        let (scalars, bases) = (scalars.clone(), bases.clone());
        c.bench_function("naive mul then add", move |b| {
            b.iter(|| apk_proofs::utils::mul_then_add(black_box(&bases), black_box(&scalars)))
        });
    }

    let nu: G::ScalarField = u128::rand(rng).into();

    {
        let bases = bases.clone();
        c.bench_function("128-bit Horner", move |b| {
            b.iter(|| apk_proofs::utils::horner(black_box(&bases), black_box(nu)))
        });
    }
}

fn bw6_subgroup_check(c: &mut Criterion) {
    let rng = &mut test_rng();

    let p = ark_bw6_761::G1Projective::rand(rng);
    let p_affine = p.into_affine();

    c.bench_function("subgroup check: mul by group order", move |b| {
        b.iter(|| (black_box(&p_affine)).is_in_correct_subgroup_assuming_on_curve())
    });

    c.bench_function("subgroup check: GLV", move |b| {
        b.iter(|| {
            apk_proofs::endo::subgroup_check::<ark_bw6_761::Config>(
                black_box(&p),
                black_box(OMEGA),
                black_box(U),
            )
        })
    });
}

fn amplification(c: &mut Criterion) {
    use apk_proofs::domains::Domains;

    let mut group = c.benchmark_group("amplification");

    let rng = &mut test_rng();
    let log_domain_size_range = vec![10, 16];

    for log_domain_size in log_domain_size_range {
        let n = 2u32.pow(log_domain_size) as usize;
        let domains = Domains::<Fr, BenchDomain>::new(n);

        let evals = (0..n).map(|_| Fr::rand(rng)).collect::<Vec<_>>();

        group.bench_with_input(
            BenchmarkId::new("2x", log_domain_size),
            &log_domain_size,
            |b, _| {
                b.iter(|| {
                    let poly = domains.interpolate(evals.clone());
                    domains.domain2x().fft(&poly.coeffs)
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("2x-coset", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| domains.amplify_x2(evals.clone())),
        );

        group.bench_with_input(
            BenchmarkId::new("4x", log_domain_size),
            &log_domain_size,
            |b, _| {
                b.iter(|| {
                    let poly = domains.interpolate(evals.clone());
                    domains.domain4x().fft(&poly.coeffs)
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("4x-coset", log_domain_size),
            &log_domain_size,
            |b, _| b.iter(|| domains.amplify_x4(evals.clone())),
        );
    }
}

fn verification(c: &mut Criterion) {
    use apk_proofs::{ApkConfig, Bitmask, Bls12_377Config, Prover, Verifier};

    let mut group = c.benchmark_group("verification");

    let rng = &mut test_rng();
    let log_domain_size_range = 8..=10;

    // Type aliases for clarity
    type InnerCurve = ark_bls12_377::G1Projective;
    type OuterCurve = ark_bw6_761::G1Projective;
    type TestPCS = KZG<BW6_761>;

    for log_domain_size in log_domain_size_range {
        let keyset_size = (2u32.pow(log_domain_size) - 1) as usize;
        let pks: Vec<InnerCurve> = (0..keyset_size).map(|_| InnerCurve::rand(rng)).collect();

        let keyset = Keyset::<InnerCurve, OuterCurve, BenchDomain>::new(pks).unwrap();

        let pcs_params = setup::params_for_keyset::<_, TestPCS, BenchDomain, _>(
            // Insecure source: fine for a benchmark.
            &mut setup::InsecureSetup::new(rng),
            keyset_size,
        )
        .unwrap();
        let pks_comm = keyset.commit::<TestPCS>(&pcs_params.ck()).unwrap();

        let bitmask = Bitmask::from_bits(&vec![true; keyset_size]);

        let prover = Prover::<InnerCurve, OuterCurve, TestPCS, BenchDomain>::new(
            keyset,
            &pks_comm,
            pcs_params.clone(),
            Bls12_377Config::transcript(),
        )
        .unwrap();

        let proof_basic = prover.prove_simple(bitmask.clone()).unwrap();
        let proof_packed = prover.prove_packed(bitmask.clone()).unwrap();
        let proof_counting = prover.prove_counting(bitmask.clone()).unwrap();

        let create_verifier = || {
            Verifier::<InnerCurve, OuterCurve, TestPCS, BenchDomain>::try_new(
                pcs_params.raw_vk(),
                pks_comm.clone(),
                Bls12_377Config::transcript(),
            )
            .unwrap()
        };

        group.bench_with_input(
            BenchmarkId::new("basic", log_domain_size),
            &log_domain_size,
            |b, _| {
                b.iter(|| {
                    let verifier = create_verifier();
                    assert!(verifier
                        .verify_simple(&proof_basic.1, &proof_basic.0)
                        .unwrap());
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("packed", log_domain_size),
            &log_domain_size,
            |b, _| {
                b.iter(|| {
                    let verifier = create_verifier();
                    assert!(verifier
                        .verify_packed(&proof_packed.1, &proof_packed.0)
                        .unwrap());
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("counting", log_domain_size),
            &log_domain_size,
            |b, _| {
                b.iter(|| {
                    let verifier = create_verifier();
                    assert!(verifier
                        .verify_counting(&proof_counting.1, &proof_counting.0)
                        .unwrap());
                })
            },
        );
    }

    group.finish();
}

fn fft<F: FftField, D: EvaluationDomain<F>>(c: &mut Criterion) {
    let mut group = c.benchmark_group("FFT");

    let rng = &mut test_rng();

    for logn in 10..=16 {
        let n = 2usize.pow(logn);
        let domain = D::new(n).unwrap();
        let poly = DensePolynomial::<F>::rand(n - 1, rng);
        let coeffs = poly.coeffs;

        group.throughput(Throughput::Elements(n as u64));
        group.bench_with_input(BenchmarkId::from_parameter(n), &n, |b, _n| {
            b.iter(|| domain.fft(&coeffs))
        });
    }
    group.finish();
}

fn components(c: &mut Criterion) {
    msm::<ark_bw6_761::G1Affine>(c, 6);
    barycentric_evaluation::<Fr>(c, 2u32.pow(10));
    bw6_subgroup_check(c);
    amplification(c);
}

fn primitives(c: &mut Criterion) {
    fft::<Fr, Radix2EvaluationDomain<Fr>>(c);
}

criterion_group!(benches, components, verification, primitives);
criterion_main!(benches);
