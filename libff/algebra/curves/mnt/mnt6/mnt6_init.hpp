/** @file
 *****************************************************************************

 Declaration of interfaces for initializing MNT6.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT6_INIT_HPP_
#define MNT6_INIT_HPP_

#include <libff/algebra/curves/mnt/mnt46_common.hpp>
#include <libff/algebra/curves/public_params.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g1.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g2.hpp>
#include <libff/algebra/fields/fp.hpp>
#include <libff/algebra/fields/fp3.hpp>
#include <libff/algebra/fields/fp6_2over3.hpp>

namespace libff {

#define mnt6_modulus_r mnt46_modulus_B
#define mnt6_modulus_q mnt46_modulus_A

const mp_size_t mnt6_r_bitcount = mnt46_B_bitcount;
const mp_size_t mnt6_q_bitcount = mnt46_A_bitcount;

const mp_size_t mnt6_r_limbs = mnt46_B_limbs;
const mp_size_t mnt6_q_limbs = mnt46_A_limbs;

extern bigint<mnt6_r_limbs> mnt6_modulus_r;
extern bigint<mnt6_q_limbs> mnt6_modulus_q;

typedef Fp_model<mnt6_r_limbs, mnt6_modulus_r> mnt6_Fr;
typedef Fp_model<mnt6_q_limbs, mnt6_modulus_q> mnt6_Fq;
typedef Fp3_model<mnt6_q_limbs, mnt6_modulus_q> mnt6_Fq3;
typedef Fp6_2over3_model<mnt6_q_limbs, mnt6_modulus_q> mnt6_Fq6;
typedef mnt6_Fq6 mnt6_GT;

// parameters for twisted short Weierstrass curve E'/Fq3 : y^2 = x^3 + (a * twist^2) * x + (b * twist^3)
extern mnt6_Fq3 mnt6_twist_coeff_a;
extern mnt6_Fq3 mnt6_twist_coeff_b;
extern mnt6_Fq mnt6_twist_mul_by_a_c0;
extern mnt6_Fq mnt6_twist_mul_by_a_c1;
extern mnt6_Fq mnt6_twist_mul_by_a_c2;
extern mnt6_Fq mnt6_twist_mul_by_b_c0;
extern mnt6_Fq mnt6_twist_mul_by_b_c1;
extern mnt6_Fq mnt6_twist_mul_by_b_c2;
extern mnt6_Fq mnt6_twist_mul_by_q_X;
extern mnt6_Fq mnt6_twist_mul_by_q_Y;

void init_mnt6_params();

class mnt6_swparams {
public:
  // G1 parameters
  typedef mnt6_Fq Fq;
  typedef mnt6_Fr Fr;
  static Fq coeff_a;
  static Fq coeff_b;
  static Fq G1_zero_X;
  static Fq G1_zero_Y;
  static Fq G1_zero_Z;
  static Fq G1_one_X;
  static Fq G1_one_Y;
  static Fq G1_one_Z;

  // G2 parameters
  typedef mnt6_Fq3 twist_field;
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
  typedef mnt6_Fq6 Fqk;
  typedef mnt6_Fq6 GT;
  static bigint<mnt6_q_limbs> ate_loop_count;
  static bool ate_is_loop_count_neg;
  static bigint<6*mnt6_q_limbs> final_exponent;
  static bigint<mnt6_q_limbs> final_exponent_last_chunk_abs_of_w0;
  static bool final_exponent_last_chunk_is_w0_neg;
  static bigint<mnt6_q_limbs> final_exponent_last_chunk_w1;

  static twist_field embed_as_first_coordinate(const Fq &x);
  static Fqk final_exponentiation_first_chunk(const Fqk &elt, const Fqk &elt_inv);
  static Fqk special_mul(const Fqk &x, const Fqk &y);
};

typedef short_weierstrass_G1<mnt6_swparams> mnt6_G1;
typedef short_weierstrass_G2<mnt6_swparams> mnt6_G2;

} // libff

#endif // MNT6_INIT_HPP_
