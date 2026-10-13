use ark_ec::CurveGroup;
use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use w3f_pcs::pcs::PCS;

use crate::domain::DomainSet;
use crate::domains::Domains;
use crate::piop::affine_addition::{
    AffineAdditionEvaluations, AffineAdditionRegisters, PartialSumsPolynomials,
};
use crate::piop::{ProverProtocol, RegisterEvaluations};
use crate::{utils, AccountablePublicInput, AccumulatorSeed, Bitmask, Keyset};

#[derive(CanonicalSerialize, CanonicalDeserialize)]
pub struct AffineAdditionEvaluationsWithoutBitmask<F: FftField> {
    pub keyset: (F, F),
    pub partial_sums: (F, F),
}

impl<F: FftField> RegisterEvaluations<F> for AffineAdditionEvaluationsWithoutBitmask<F> {
    fn as_vec(&self) -> Vec<F> {
        vec![
            self.keyset.0,
            self.keyset.1,
            self.partial_sums.0,
            self.partial_sums.1,
        ]
    }
}

pub struct BasicRegisterBuilder<F: PrimeField, D: DomainSet<F>> {
    registers: AffineAdditionRegisters<F, D>,
    register_evaluations: Option<AffineAdditionEvaluations<F>>,
}

impl<IC, OC, S, D> ProverProtocol<IC, OC, S, D> for BasicRegisterBuilder<OC::ScalarField, D>
where
    IC: AccumulatorSeed,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField>,
    S: PCS<OC::ScalarField>,
    D: DomainSet<OC::ScalarField>,
{
    type P1 = PartialSumsPolynomials<OC::ScalarField>;
    type P2 = ();
    type E = AffineAdditionEvaluationsWithoutBitmask<OC::ScalarField>;
    type PI = AccountablePublicInput<IC>;

    fn init(
        domains: Domains<OC::ScalarField, D>,
        bitmask: Bitmask,
        keyset: Keyset<IC, OC, D>,
    ) -> Self {
        BasicRegisterBuilder {
            registers: AffineAdditionRegisters::<OC::ScalarField, D>::new(
                domains,
                keyset,
                &bitmask.to_bits(),
            ),
            register_evaluations: None,
        }
    }

    fn get_register_polynomials_to_commit1(&self) -> PartialSumsPolynomials<OC::ScalarField> {
        let polys = self.registers.get_register_polynomials();
        polys.partial_sums
    }

    fn get_register_polynomials_to_commit2(&mut self, _verifier_challenge: OC::ScalarField) -> () {
        ()
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<OC::ScalarField>> {
        let polys = self.registers.get_register_polynomials();
        [polys.keyset, polys.partial_sums].concat()
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<OC::ScalarField>> {
        self.registers.compute_constraint_polynomials::<IC, OC>()
    }

    // bitmask register polynomial is not committed to...
    fn evaluate_register_polynomials(
        &mut self,
        point: OC::ScalarField,
    ) -> AffineAdditionEvaluationsWithoutBitmask<OC::ScalarField> {
        let evals: AffineAdditionEvaluations<OC::ScalarField> =
            self.registers.evaluate_register_polynomials(point);
        self.register_evaluations = Some(evals.clone());
        AffineAdditionEvaluationsWithoutBitmask {
            keyset: evals.keyset,
            partial_sums: evals.partial_sums,
        }
    }

    fn compute_linearization_polynomial(
        &self,
        phi: OC::ScalarField,
        zeta: OC::ScalarField,
    ) -> DensePolynomial<OC::ScalarField> {
        #[allow(clippy::expect_used, reason = "see ROUNDS_IN_ORDER")]
        let evals = self
            .register_evaluations
            .as_ref()
            .expect(crate::piop::ROUNDS_IN_ORDER);
        let parts = self.registers.compute_constraints_linearized(evals, zeta);
        utils::randomize(phi, &parts)
    }
}
