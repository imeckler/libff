/** @file
 *****************************************************************************

 Declaration of interfaces for initializing MNT4.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT4_INIT_HPP_
#define MNT4_INIT_HPP_

#include <libff/algebra/curves/mnt/mnt46_common.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g1.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g2.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp2.hpp>
#include <libff/algebra/fields/fp4.hpp>

namespace libff {

#define mnt4_modulus_r mnt46_modulus_A
#define mnt4_modulus_q mnt46_modulus_B

const mp_size_t mnt4_r_bitcount = mnt46_A_bitcount;
const mp_size_t mnt4_q_bitcount = mnt46_B_bitcount;

const mp_size_t mnt4_r_limbs = mnt46_A_limbs;
const mp_size_t mnt4_q_limbs = mnt46_B_limbs;

extern bigint<mnt4_r_limbs> mnt4_modulus_r;
extern bigint<mnt4_q_limbs> mnt4_modulus_q;

typedef Fp_model<mnt4_r_limbs, mnt4_modulus_r> mnt4_Fr;
typedef Fp_model<mnt4_q_limbs, mnt4_modulus_q> mnt4_Fq;
typedef Fp2_model<mnt4_q_limbs, mnt4_modulus_q> mnt4_Fq2;
typedef Fp4_model<mnt4_q_limbs, mnt4_modulus_q> mnt4_Fq4;
typedef mnt4_Fq4 mnt4_GT;

// parameters for twisted short Weierstrass curve E'/Fq2 : y^2 = x^3 + (a * twist^2) * x + (b * twist^3)
// TODO: Delete extern mnt4_Fq2 mnt4_twist;
extern mnt4_Fq2 mnt4_twist_coeff_a;
extern mnt4_Fq2 mnt4_twist_coeff_b;
extern mnt4_Fq mnt4_twist_mul_by_a_c0;
extern mnt4_Fq mnt4_twist_mul_by_a_c1;
extern mnt4_Fq mnt4_twist_mul_by_b_c0;
extern mnt4_Fq mnt4_twist_mul_by_b_c1;
extern mnt4_Fq mnt4_twist_mul_by_q_X;
extern mnt4_Fq mnt4_twist_mul_by_q_Y;

void init_mnt4_params();

class mnt4_swparams {
public:
  // G1 parameters
  typedef mnt4_Fq Fq;
  typedef mnt4_Fr Fr;
  static Fq coeff_a;
  static Fq coeff_b;
  static Fq G1_zero_X;
  static Fq G1_zero_Y;
  static Fq G1_zero_Z;
  static Fq G1_one_X;
  static Fq G1_one_Y;
  static Fq G1_one_Z;

  // G2 parameters
  typedef mnt4_Fq2 twist_field;
  static twist_field twist;

  static twist_field twist_coeff_a;
  static twist_field twist_coeff_b;

  static twist_field G2_zero_X;
  static twist_field G2_zero_Y;
  static twist_field G2_zero_Z;
  static twist_field G2_one_X;
  static twist_field G2_one_Y;
  static twist_field G2_one_Z;

  static twist_field mul_by_a(const twist_field &elt);
  static twist_field mul_by_b(const twist_field &elt);

  static Fq twist_mul_by_q_X;
  static Fq twist_mul_by_q_Y;

  // pairing parameters
  typedef mnt4_Fq4 Fqk;
  typedef mnt4_Fq4 GT;
  static bigint<mnt4_q_limbs> ate_loop_count;
  static bool ate_is_loop_count_neg;
  static bigint<4*mnt4_q_limbs> final_exponent;
  static bigint<mnt4_q_limbs> final_exponent_last_chunk_abs_of_w0;
  static bool final_exponent_last_chunk_is_w0_neg;
  static bigint<mnt4_q_limbs> final_exponent_last_chunk_w1;

  static Fqk final_exponentiation_first_chunk(const Fqk &elt, const Fqk &elt_inv);
  static Fqk special_mul(const Fqk &x, const Fqk &y);
};

typedef short_weierstrass_G1<mnt4_swparams> mnt4_G1;
typedef short_weierstrass_G2<mnt4_swparams> mnt4_G2;

} // libff

#endif // MNT4_INIT_HPP_
