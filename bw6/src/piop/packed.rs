use ark_ec::CurveGroup;
use ark_ff::{FftField, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use w3f_pcs::pcs::PCS;

use crate::domain::DomainSet;
use crate::domains::Domains;
use crate::piop::affine_addition::{AffineAdditionRegisters, PartialSumsAndBitmaskPolynomials};
use crate::piop::bitmask_packing::{
    BitmaskPackingPolynomials, BitmaskPackingRegisters, SuccinctAccountableRegisterEvaluations,
};
use crate::piop::ProverProtocol;
use crate::{utils, AccountablePublicInput, AccumulatorSeed, Bitmask, Keyset};

pub struct PackedRegisterBuilder<F: PrimeField, D: DomainSet<F>> {
    bitmask: Bitmask,
    affine_addition_registers: AffineAdditionRegisters<F, D>,
    bitmask_packing_registers: Option<BitmaskPackingRegisters<F, D>>,
    register_evaluations: Option<SuccinctAccountableRegisterEvaluations<F>>,
}

#[allow(clippy::expect_used, reason = "see ROUNDS_IN_ORDER")]
impl<IC, OC, S, D> ProverProtocol<IC, OC, S, D> for PackedRegisterBuilder<OC::ScalarField, D>
where
    IC: AccumulatorSeed,
    OC: CurveGroup,
    OC::ScalarField: From<IC::BaseField> + FftField,
    S: PCS<OC::ScalarField>,
    D: DomainSet<OC::ScalarField>,
{
    type P1 = PartialSumsAndBitmaskPolynomials<OC::ScalarField>;
    type P2 = BitmaskPackingPolynomials<OC::ScalarField>;
    type E = SuccinctAccountableRegisterEvaluations<OC::ScalarField>;
    type PI = AccountablePublicInput<IC>;

    fn init(
        domains: Domains<OC::ScalarField, D>,
        bitmask: Bitmask,
        keyset: Keyset<IC, OC, D>,
    ) -> Self {
        PackedRegisterBuilder {
            bitmask: bitmask.clone(),
            affine_addition_registers: AffineAdditionRegisters::new(
                domains,
                keyset,
                &bitmask.to_bits(),
            ),
            bitmask_packing_registers: None,
            register_evaluations: None,
        }
    }

    fn get_register_polynomials_to_commit1(
        &self,
    ) -> PartialSumsAndBitmaskPolynomials<OC::ScalarField> {
        self.affine_addition_registers
            .get_partial_sums_and_bitmask_polynomials()
    }

    fn get_register_polynomials_to_commit2(
        &mut self,
        bitmask_chunks_aggregation_challenge: OC::ScalarField,
    ) -> BitmaskPackingPolynomials<OC::ScalarField> {
        let bitmask_packing_registers = BitmaskPackingRegisters::new(
            self.affine_addition_registers.domains.clone(),
            &self.bitmask,
            bitmask_chunks_aggregation_challenge,
        );
        let res = bitmask_packing_registers.get_register_polynomials();
        self.bitmask_packing_registers = Some(bitmask_packing_registers);
        res
    }

    fn get_register_polynomials_to_open(self) -> Vec<DensePolynomial<OC::ScalarField>> {
        let affine_addition_polys = self
            .affine_addition_registers
            .get_register_polynomials()
            .to_vec();
        let bitmask_packing_polys = self
            .bitmask_packing_registers
            .expect(crate::piop::ROUNDS_IN_ORDER)
            .get_register_polynomials()
            .to_vec();
        let mut polys = vec![];
        polys.extend(affine_addition_polys);
        polys.extend(bitmask_packing_polys);
        polys
    }

    fn compute_constraint_polynomials(&self) -> Vec<DensePolynomial<OC::ScalarField>> {
        let affine_addition_constraints = self
            .affine_addition_registers
            .compute_constraint_polynomials::<IC, OC>();
        let bitmask_packing_constraints = self
            .bitmask_packing_registers
            .as_ref()
            .expect(crate::piop::ROUNDS_IN_ORDER)
            .compute_constraint_polynomials();
        let mut constraints = vec![];
        constraints.extend(affine_addition_constraints);
        constraints.extend(bitmask_packing_constraints);
        constraints
    }

    fn evaluate_register_polynomials(
        &mut self,
        point: OC::ScalarField,
    ) -> SuccinctAccountableRegisterEvaluations<OC::ScalarField> {
        let affine_addition_evals = self
            .affine_addition_registers
            .evaluate_register_polynomials(point);
        let bitmask_packing_evals = self
            .bitmask_packing_registers
            .as_ref()
            .expect(crate::piop::ROUNDS_IN_ORDER)
            .evaluate_register_polynomials(point);
        let evals = SuccinctAccountableRegisterEvaluations {
            c: bitmask_packing_evals.0,
            acc: bitmask_packing_evals.1,
            basic_evaluations: affine_addition_evals,
        };
        self.register_evaluations = Some(evals.clone());
        evals
    }

    fn compute_linearization_polynomial(
        &self,
        phi: OC::ScalarField,
        zeta: OC::ScalarField,
    ) -> DensePolynomial<OC::ScalarField> {
        let evals = self
            .register_evaluations
            .as_ref()
            .expect(crate::piop::ROUNDS_IN_ORDER);

        let affine_addition_parts = self
            .affine_addition_registers
            .compute_constraints_linearized(&evals.basic_evaluations, zeta);
        let bitmask_packing_parts = self
            .bitmask_packing_registers
            .as_ref()
            .expect(crate::piop::ROUNDS_IN_ORDER)
            .compute_constraints_linearized();

        let mut parts = vec![];
        parts.extend(affine_addition_parts);
        parts.extend(bitmask_packing_parts);
        utils::randomize(phi, &parts)
    }
}
