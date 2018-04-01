/** @file
 *****************************************************************************

 Declaration of interfaces for pairing operations on MNT4.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT4_PAIRING_HPP_
#define MNT4_PAIRING_HPP_

#include <vector>

#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_pairing.hpp>

namespace libff {

typedef short_weierstrass_G1_precomp<mnt4_swparams> mnt4_G1_precomp;
typedef short_weierstrass_G2_precomp<mnt4_swparams> mnt4_G2_precomp;

typedef short_weierstrass_affine_ate_G1_precomputation<mnt4_swparams> mnt4_affine_ate_G1_precomputation;
typedef short_weierstrass_affine_ate_G2_precomputation<mnt4_swparams> mnt4_affine_ate_G2_precomputation;

typedef short_weierstrass_affine_ate_G1_precomputation<mnt4_swparams> mnt4_affine_ate_G1_precomputation;
typedef short_weierstrass_affine_ate_coeffs<mnt4_swparams> mnt4_affine_ate_coeffs;
typedef short_weierstrass_affine_ate_G2_precomputation<mnt4_swparams> mnt4_affine_ate_G2_precomputation;

typedef short_weierstrass_ate_G1_precomp<mnt4_swparams> mnt4_ate_G1_precomp;
typedef short_weierstrass_ate_G2_precomp<mnt4_swparams> mnt4_ate_G2_precomp;
typedef short_weierstrass_ate_dbl_coeffs<mnt4_swparams> mnt4_ate_dbl_coeffs;
typedef short_weierstrass_ate_add_coeffs<mnt4_swparams> mnt4_ate_add_coeffs;

/* final exponentiation */

mnt4_Fq4 mnt4_final_exponentiation_last_chunk(const mnt4_Fq4 &elt,
                                              const mnt4_Fq4 &elt_inv);
mnt4_Fq4 mnt4_final_exponentiation_first_chunk(const mnt4_Fq4 &elt,
                                               const mnt4_Fq4 &elt_inv);
mnt4_GT mnt4_final_exponentiation(const mnt4_Fq4 &elt);

/* affine ate miller loop */
mnt4_affine_ate_G1_precomputation mnt4_affine_ate_precompute_G1(const mnt4_G1& P);
mnt4_affine_ate_G2_precomputation mnt4_affine_ate_precompute_G2(const mnt4_G2& Q);

mnt4_Fq4 mnt4_affine_ate_miller_loop(const mnt4_affine_ate_G1_precomputation &prec_P,
                                     const mnt4_affine_ate_G2_precomputation &prec_Q);

/* ate pairing */
mnt4_ate_G1_precomp mnt4_ate_precompute_G1(const mnt4_G1& P);
mnt4_ate_G2_precomp mnt4_ate_precompute_G2(const mnt4_G2& Q);

mnt4_Fq4 mnt4_ate_miller_loop(const mnt4_ate_G1_precomp &prec_P,
                                    const mnt4_ate_G2_precomp &prec_Q);
mnt4_Fq4 mnt4_ate_double_miller_loop(const mnt4_ate_G1_precomp &prec_P1,
                                           const mnt4_ate_G2_precomp &prec_Q1,
                                           const mnt4_ate_G1_precomp &prec_P2,
                                           const mnt4_ate_G2_precomp &prec_Q2);

mnt4_Fq4 mnt4_ate_pairing(const mnt4_G1& P,
                          const mnt4_G2 &Q);
mnt4_GT mnt4_ate_reduced_pairing(const mnt4_G1 &P,
                                 const mnt4_G2 &Q);

/* choice of pairing */
typedef mnt4_ate_G1_precomp mnt4_G1_precomp;
typedef mnt4_ate_G2_precomp mnt4_G2_precomp;

mnt4_G1_precomp mnt4_precompute_G1(const mnt4_G1& P);

mnt4_G2_precomp mnt4_precompute_G2(const mnt4_G2& Q);

mnt4_Fq4 mnt4_miller_loop(const mnt4_G1_precomp &prec_P,
                          const mnt4_G2_precomp &prec_Q);

mnt4_Fq4 mnt4_double_miller_loop(const mnt4_G1_precomp &prec_P1,
                                 const mnt4_G2_precomp &prec_Q1,
                                 const mnt4_G1_precomp &prec_P2,
                                 const mnt4_G2_precomp &prec_Q2);

mnt4_Fq4 mnt4_pairing(const mnt4_G1& P,
                      const mnt4_G2 &Q);

mnt4_GT mnt4_reduced_pairing(const mnt4_G1 &P,
                             const mnt4_G2 &Q);

mnt4_GT mnt4_affine_reduced_pairing(const mnt4_G1 &P,
                                    const mnt4_G2 &Q);


} // libff

#endif // MNT4_PAIRING_HPP_
