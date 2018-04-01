/** @file
 *****************************************************************************

 Implementation of interfaces for pairing operations on MNT4.

 See mnt4_pairing.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pairing.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

/* final exponentiation */
mnt4_Fq4 mnt4_final_exponentiation_last_chunk(const mnt4_Fq4 &elt,
                                              const mnt4_Fq4 &elt_inv)
{
    return short_weierstrass_final_exponentiation_last_chunk<mnt4_swparams>(elt, elt_inv);
}

mnt4_Fq4 mnt4_final_exponentiation_first_chunk(const mnt4_Fq4 &elt,
                                               const mnt4_Fq4 &elt_inv)
{
    return short_weierstrass_final_exponentiation_first_chunk<mnt4_swparams>(elt, elt_inv);
}

mnt4_GT mnt4_final_exponentiation(const mnt4_Fq4 &elt)
{
    return short_weierstrass_final_exponentiation<mnt4_swparams>(elt);
}

/* affine ate miller loop */
mnt4_affine_ate_G1_precomputation mnt4_affine_ate_precompute_G1(const mnt4_G1& P)
{
    return short_weierstrass_affine_ate_precompute_G1<mnt4_swparams>(P);
}
mnt4_affine_ate_G2_precomputation mnt4_affine_ate_precompute_G2(const mnt4_G2& Q)
{
    return short_weierstrass_affine_ate_precompute_G2<mnt4_swparams>(Q);
}

mnt4_Fq4 mnt4_affine_ate_miller_loop(const mnt4_affine_ate_G1_precomputation &prec_P,
                                     const mnt4_affine_ate_G2_precomputation &prec_Q)
{
    return short_weierstrass_affine_ate_miller_loop<mnt4_swparams>(prec_P, prec_Q);
}

/* ate pairing */
mnt4_ate_G1_precomp mnt4_ate_precompute_G1(const mnt4_G1& P)
{
    return short_weierstrass_ate_precompute_G1<mnt4_swparams>(P);
}
mnt4_ate_G2_precomp mnt4_ate_precompute_G2(const mnt4_G2& Q)
{
    return short_weierstrass_ate_precompute_G2<mnt4_swparams>(Q);
}

mnt4_Fq4 mnt4_ate_miller_loop(const mnt4_ate_G1_precomp &prec_P,
                                    const mnt4_ate_G2_precomp &prec_Q)
{
    return short_weierstrass_ate_miller_loop<mnt4_swparams>(prec_P, prec_Q);
}
mnt4_Fq4 mnt4_ate_double_miller_loop(const mnt4_ate_G1_precomp &prec_P1,
                                           const mnt4_ate_G2_precomp &prec_Q1,
                                           const mnt4_ate_G1_precomp &prec_P2,
                                           const mnt4_ate_G2_precomp &prec_Q2)
{
    return short_weierstrass_ate_double_miller_loop<mnt4_swparams>(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

mnt4_Fq4 mnt4_ate_pairing(const mnt4_G1& P,
                          const mnt4_G2 &Q)
{
}
mnt4_GT mnt4_ate_reduced_pairing(const mnt4_G1 &P,
                                 const mnt4_G2 &Q)
{
    return short_weierstrass_ate_reduced_pairing<mnt4_swparams>(P, Q);
}

/* choice of pairing */

mnt4_G1_precomp mnt4_precompute_G1(const mnt4_G1& P)
{
    return short_weierstrass_precompute_G1<mnt4_swparams>(P);
}

mnt4_G2_precomp mnt4_precompute_G2(const mnt4_G2& Q)
{
    return short_weierstrass_precompute_G2<mnt4_swparams>(Q);
}

mnt4_Fq4 mnt4_miller_loop(const mnt4_G1_precomp &prec_P,
                          const mnt4_G2_precomp &prec_Q)
{
    return short_weierstrass_miller_loop<mnt4_swparams>(prec_P, prec_Q);
}

mnt4_Fq4 mnt4_double_miller_loop(const mnt4_G1_precomp &prec_P1,
                                 const mnt4_G2_precomp &prec_Q1,
                                 const mnt4_G1_precomp &prec_P2,
                                 const mnt4_G2_precomp &prec_Q2)
{
    return short_weierstrass_double_miller_loop<mnt4_swparams>(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

mnt4_Fq4 mnt4_pairing(const mnt4_G1& P,
                      const mnt4_G2 &Q)
{
    return short_weierstrass_pairing<mnt4_swparams>(P, Q);
}

mnt4_GT mnt4_reduced_pairing(const mnt4_G1 &P,
                             const mnt4_G2 &Q)
{
    return short_weierstrass_reduced_pairing<mnt4_swparams>(P, Q);
}

mnt4_GT mnt4_affine_reduced_pairing(const mnt4_G1 &P,
                                    const mnt4_G2 &Q)
{
    return short_weierstrass_affine_reduced_pairing<mnt4_swparams>(P, Q);
}

} // libff
