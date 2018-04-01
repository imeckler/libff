/** @file
 *****************************************************************************

 Declaration of interfaces for pairing operations on MNT6.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT6_PAIRING_HPP_
#define MNT6_PAIRING_HPP_

#include <vector>

#include <libff/algebra/curves/mnt/mnt6/mnt6_init.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_pairing.hpp>

namespace libff {

/* final exponentiation */

mnt6_Fq6 mnt6_final_exponentiation_last_chunk(const mnt6_Fq6 &elt,
                                              const mnt6_Fq6 &elt_inv);
mnt6_Fq6 mnt6_final_exponentiation_first_chunk(const mnt6_Fq6 &elt,
                                               const mnt6_Fq6 &elt_inv);
mnt6_GT mnt6_final_exponentiation(const mnt6_Fq6 &elt);

/* affine ate miller loop */

typedef short_weierstrass_affine_ate_G1_precomputation<mnt6_swparams> mnt6_affine_ate_G1_precomputation;
typedef short_weierstrass_affine_ate_coeffs<mnt6_swparams> mnt6_affine_ate_coeffs;
typedef short_weierstrass_affine_ate_G2_precomputation<mnt6_swparams> mnt6_affine_ate_G2_precomputation;

mnt6_affine_ate_G1_precomputation mnt6_affine_ate_precompute_G1(const mnt6_G1& P);
mnt6_affine_ate_G2_precomputation mnt6_affine_ate_precompute_G2(const mnt6_G2& Q);

mnt6_Fq6 mnt6_affine_ate_miller_loop(const mnt6_affine_ate_G1_precomputation &prec_P,
                                     const mnt6_affine_ate_G2_precomputation &prec_Q);

/* ate pairing */

typedef short_weierstrass_ate_G1_precomp<mnt6_swparams> mnt6_ate_G1_precomp;
typedef short_weierstrass_ate_dbl_coeffs<mnt6_swparams> mnt6_ate_dbl_coeffs;
typedef short_weierstrass_ate_add_coeffs<mnt6_swparams> mnt6_ate_add_coeffs;
typedef short_weierstrass_ate_G2_precomp<mnt6_swparams> mnt6_ate_G2_precomp;

mnt6_ate_G1_precomp mnt6_ate_precompute_G1(const mnt6_G1& P);
mnt6_ate_G2_precomp mnt6_ate_precompute_G2(const mnt6_G2& Q);

mnt6_Fq6 mnt6_ate_miller_loop(const mnt6_ate_G1_precomp &prec_P,
                              const mnt6_ate_G2_precomp &prec_Q);
mnt6_Fq6 mnt6_ate_double_miller_loop(const mnt6_ate_G1_precomp &prec_P1,
                                     const mnt6_ate_G2_precomp &prec_Q1,
                                     const mnt6_ate_G1_precomp &prec_P2,
                                     const mnt6_ate_G2_precomp &prec_Q2);

mnt6_Fq6 mnt6_ate_pairing(const mnt6_G1& P,
                          const mnt6_G2 &Q);
mnt6_GT mnt6_ate_reduced_pairing(const mnt6_G1 &P,
                                 const mnt6_G2 &Q);

/* choice of pairing */

typedef mnt6_ate_G1_precomp mnt6_G1_precomp;
typedef mnt6_ate_G2_precomp mnt6_G2_precomp;

mnt6_G1_precomp mnt6_precompute_G1(const mnt6_G1& P);

mnt6_G2_precomp mnt6_precompute_G2(const mnt6_G2& Q);

mnt6_Fq6 mnt6_miller_loop(const mnt6_G1_precomp &prec_P,
                          const mnt6_G2_precomp &prec_Q);

mnt6_Fq6 mnt6_double_miller_loop(const mnt6_G1_precomp &prec_P1,
                                 const mnt6_G2_precomp &prec_Q1,
                                 const mnt6_G1_precomp &prec_P2,
                                 const mnt6_G2_precomp &prec_Q2);

mnt6_Fq6 mnt6_pairing(const mnt6_G1& P,
                      const mnt6_G2 &Q);

mnt6_GT mnt6_reduced_pairing(const mnt6_G1 &P,
                             const mnt6_G2 &Q);

mnt6_GT mnt6_affine_reduced_pairing(const mnt6_G1 &P,
                                    const mnt6_G2 &Q);

} // libff

#endif // MNT6_PAIRING_HPP_
