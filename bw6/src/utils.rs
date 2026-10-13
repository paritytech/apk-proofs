use ark_ec::{AffineRepr, CurveGroup};
use ark_ff::{batch_inversion, FftField, PrimeField};
use ark_ff::{Field, Zero};
use ark_poly::Polynomial;

use crate::domain::FftDomain;
use crate::Bitmask;

// Evaluates a polynomial represented as evaluations over a multiplicative subgroup (aka in Lagrange basis) at a point.
// f = sum(fi * Li), where Li is the i-th Lagrange basis polynomial, and fi = f(w^i)
// Li(z) = (z^n -1)/n * li(z), where li(z) := w^i /(z-w^i), see https://hackmd.io/xTta-c--SFyOv9Kl3Q9Jjw
// Then f(z) = sum(fi * Li(z)) = (z^n -1)/n sum(fi * li(z))
// To avoid inversions when computing li(z) we instead compute li_inv = (z-w^i)/w^i = z/w^i - 1 = z * w_inv^i - 1
// using n multiplications to accumulate z * w_inv^i and then perform batch inversion.
// Batch inversion costs 1 inv + 3n muls, and n more muls is required for fi * li(z),
// resulting in around 1 inv + 5n muls
pub fn barycentric_eval_at<F: PrimeField, D: FftDomain<F>>(z: F, evals: &Vec<F>, domain: &D) -> F {
    let n = domain.size();
    assert_eq!(evals.len(), n);
    // Generic exponentiation rather than a chain of squarings: n is not a power of two in
    // general. For radix-2 domains this computes exactly the same value.
    let mut z_n = z.pow([n as u64]);
    z_n -= F::one();
    z_n *= &domain.size_inv(); // (z^n-1)/n

    let mut li_inv = Vec::with_capacity(n);
    let mut acc = z;
    for _ in 0..n {
        li_inv.push(acc - F::one());
        acc *= domain.generator_inv();
    }
    batch_inversion(&mut li_inv);
    let s = evals.iter().zip(li_inv).map(|(fi, li)| li * fi).sum::<F>();
    z_n * s
}

pub fn barycentric_eval_binary_at<F: PrimeField, D: FftDomain<F>>(
    z: F,
    evals: &Bitmask,
    domain: &D,
) -> F {
    let mut z_n = z.pow([domain.size() as u64]);
    z_n -= F::one();
    z_n *= &domain.size_inv(); // (z^n-1)/n

    let mut li_inv = Vec::with_capacity(evals.count_ones());
    let mut acc = z;
    for b in evals.to_bits() {
        if b {
            li_inv.push(acc - F::one());
        }
        acc *= domain.generator_inv();
    }

    batch_inversion(&mut li_inv);
    let s: F = li_inv.iter().sum();
    z_n * s
}

/// Values of the polynomials at a point z
pub struct LagrangeEvaluations<F: FftField> {
    pub vanishing_polynomial: F, // z^n - 1
    pub l_first: F,              // L_0(z)
    pub l_last: F,               // L_{n-1}(z)
    pub zeta_minus_omega_inv: F, // z - \omega^{-1}
    pub zeta_omega: F,           // z * \omega
}

//TODO: move to domains
pub fn lagrange_evaluations<F: PrimeField, D: FftDomain<F>>(
    z: F,
    domain: &D,
) -> LagrangeEvaluations<F> {
    // TODO: reuse this code with barycentric_eval methods
    let z_n = z.pow([domain.size() as u64]);

    let z_n_minus_one = z_n - F::one();
    let z_n_minus_one_div_n = z_n_minus_one * domain.size_inv();

    let mut inv = [z - F::one(), domain.generator() * z - F::one()];
    batch_inversion(&mut inv);
    LagrangeEvaluations {
        vanishing_polynomial: z_n_minus_one,
        l_first: z_n_minus_one_div_n * inv[0],
        l_last: z_n_minus_one_div_n * inv[1],
        zeta_minus_omega_inv: z - domain.generator_inv(),
        zeta_omega: z * domain.generator(),
    }
}

pub fn mul_then_add<G: AffineRepr>(bases: &[G], scalars: &[G::ScalarField]) -> G::Group {
    bases.iter().zip(scalars).map(|(&b, s)| b * s).sum()
}

pub fn horner<G: AffineRepr>(bases: &[G], nu: G::ScalarField) -> G {
    bases
        .iter()
        .rev()
        .fold(G::Group::zero(), |acc, b| (acc * nu) + b)
        .into_affine()
}

pub fn horner_field<F: Field>(bases: &[F], nu: F) -> F {
    bases.iter().rev().fold(F::zero(), |acc, b| nu * acc + b)
}

/// (max_exp+1)-sized vec: 1, base, base^2,... ,base^{max_exp}
pub fn powers<F: Field>(base: F, max_exp: usize) -> Vec<F> {
    let mut result = Vec::with_capacity(max_exp + 1);
    result.push(F::one());
    if max_exp > 0 {
        result.push(base);
    }
    let mut curr = base;
    for _ in 1..max_exp {
        curr *= base;
        result.push(curr);
    }
    result
}

pub fn randomize<P, F>(r: F, polys: &[P]) -> P
where
    F: Field,
    P: Polynomial<F>,
{
    let mut res = P::zero();
    if polys.is_empty() {
        return res;
    }
    let powers = powers(r, polys.len() - 1);

    powers.into_iter().zip(polys).for_each(|(r, p)| {
        res += (r, p);
    });
    res
}

#[cfg(test)]
mod tests {
    use ark_ff::{Field, One};
    use ark_poly::{EvaluationDomain, Evaluations, Polynomial, Radix2EvaluationDomain};
    use ark_std::convert::TryInto;
    use ark_std::{test_rng, UniformRand};

    use crate::test_helpers::_random_bits;

    use super::*;

    #[test]
    pub fn test_barycentric_eval() {
        let rng = &mut test_rng();
        let n = 2u32.pow(16);
        let domain =
            Radix2EvaluationDomain::new(std::convert::TryInto::try_into(n).unwrap()).unwrap();
        let z = ark_bw6_761::Fr::rand(rng);

        let evals = (0..n)
            .map(|_| ark_bw6_761::Fr::rand(rng))
            .collect::<Vec<_>>();
        let poly = Evaluations::from_vec_and_domain(evals.clone(), domain).interpolate();
        let poly_at_z = poly.evaluate(&z);
        let rdomain = crate::Radix2Domain(domain);
        assert_eq!(barycentric_eval_at(z, &evals, &rdomain), poly_at_z);

        let bitmask = Bitmask::from_bits(&_random_bits(n.try_into().unwrap(), 1.0 / 2.0, rng));
        let bits_as_field_elements = bitmask
            .to_bits()
            .iter()
            .map(|b| {
                if *b {
                    ark_bw6_761::Fr::one()
                } else {
                    ark_bw6_761::Fr::zero()
                }
            })
            .collect::<Vec<_>>();
        let bits_poly =
            Evaluations::from_vec_and_domain(bits_as_field_elements.clone(), domain).interpolate();
        let bits_poly_at_z = bits_poly.evaluate(&z);
        assert_eq!(
            barycentric_eval_at(z, &bits_as_field_elements, &rdomain),
            bits_poly_at_z
        );
        assert_eq!(
            barycentric_eval_binary_at(z, &bitmask, &rdomain),
            bits_poly_at_z
        );
    }

    #[test]
    pub fn test_horner() {
        let n = 10;

        let rng = &mut test_rng();

        let nu = ark_bw6_761::Fr::rand(rng);
        let bases = (0..n)
            .map(|_| ark_bw6_761::G1Projective::rand(rng).into_affine())
            .collect::<Vec<_>>();

        let powers = (0..n).map(|i| nu.pow([i as u64])).collect::<Vec<_>>();

        assert_eq!(horner(&bases, nu), mul_then_add(&bases, &powers));
    }

    #[test]
    fn test_lagrange_evaluations() {
        let rng = &mut test_rng();

        let n = 16;
        let domain = Radix2EvaluationDomain::<ark_bw6_761::Fr>::new(n).unwrap();

        let z = ark_bw6_761::Fr::rand(rng);
        let evals = lagrange_evaluations(z, &crate::Radix2Domain(domain));
        assert_eq!(
            evals.vanishing_polynomial,
            domain.evaluate_vanishing_polynomial(z)
        );
        let coeffs = domain.evaluate_all_lagrange_coefficients(z);
        assert_eq!(evals.l_first, coeffs[0]);
        assert_eq!(evals.l_last, coeffs[n - 1]);
    }

    #[test]
    fn test_powers() {
        let rng = &mut test_rng();
        let base = ark_bw6_761::Fr::rand(rng);

        assert_eq!(powers(base, 0), vec![ark_bw6_761::Fr::one()]);
        assert_eq!(powers(base, 1), vec![ark_bw6_761::Fr::one(), base]);

        let max_exp = 255;
        let res = powers(base, max_exp);
        assert_eq!(res.len(), max_exp + 1);
        assert_eq!(res[0], ark_bw6_761::Fr::one());
        assert_eq!(res[1], base);
        assert_eq!(res[255], base.pow([255u64]));
    }
}
