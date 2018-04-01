/** @file
 *****************************************************************************

 Implementation of interfaces for pairing operations on MNT6.

 See mnt6_pairing.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <cassert>

#include <libff/algebra/curves/mnt/mnt6/mnt6_init.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pairing.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

/* final exponentiation */

mnt6_Fq6 mnt6_final_exponentiation_last_chunk(const mnt6_Fq6 &elt,
                                              const mnt6_Fq6 &elt_inv)
{
    return short_weierstrass_final_exponentiation_last_chunk<mnt6_swparams>(elt, elt_inv);
}

mnt6_Fq6 mnt6_final_exponentiation_first_chunk(const mnt6_Fq6 &elt,
                                               const mnt6_Fq6 &elt_inv)
{
    return short_weierstrass_final_exponentiation_first_chunk<mnt6_swparams>(elt, elt_inv);
}

mnt6_GT mnt6_final_exponentiation(const mnt6_Fq6 &elt)
{
    return short_weierstrass_final_exponentiation<mnt6_swparams>(elt);
}

/* affine ate miller loop */

mnt6_affine_ate_G1_precomputation mnt6_affine_ate_precompute_G1(const mnt6_G1& P)
{
    return short_weierstrass_affine_ate_precompute_G1(P);
}

mnt6_affine_ate_G2_precomputation mnt6_affine_ate_precompute_G2(const mnt6_G2& Q)
{
    return short_weierstrass_affine_ate_precompute_G2(Q);
}


mnt6_Fq6 mnt6_affine_ate_miller_loop(const mnt6_affine_ate_G1_precomputation &prec_P,
                                     const mnt6_affine_ate_G2_precomputation &prec_Q)
{
    return short_weierstrass_affine_ate_miller_loop(prec_P, prec_Q);
}

/* ate pairing */
mnt6_ate_G1_precomp mnt6_ate_precompute_G1(const mnt6_G1& P)
{
    return short_weierstrass_ate_precompute_G1<mnt6_swparams>(P);
}

mnt6_ate_G2_precomp mnt6_ate_precompute_G2(const mnt6_G2& Q)
{
    return short_weierstrass_ate_precompute_G2<mnt6_swparams>(Q);
}

mnt6_Fq6 mnt6_ate_miller_loop(const mnt6_ate_G1_precomp &prec_P,
                              const mnt6_ate_G2_precomp &prec_Q)
{
    return short_weierstrass_ate_miller_loop<mnt6_swparams>(prec_P, prec_Q);
}

mnt6_Fq6 mnt6_ate_double_miller_loop(const mnt6_ate_G1_precomp &prec_P1,
                                     const mnt6_ate_G2_precomp &prec_Q1,
                                     const mnt6_ate_G1_precomp &prec_P2,
                                     const mnt6_ate_G2_precomp &prec_Q2)
{
    return short_weierstrass_ate_double_miller_loop<mnt6_swparams>(prec_P1, prec_Q1, prec_P2, prec_Q2);
}


mnt6_Fq6 mnt6_ate_pairing(const mnt6_G1& P,
                          const mnt6_G2 &Q)
{
    return short_weierstrass_ate_pairing<mnt6_swparams>(P, Q);
}

mnt6_GT mnt6_ate_reduced_pairing(const mnt6_G1 &P,
                                 const mnt6_G2 &Q)
{
    return short_weierstrass_ate_reduced_pairing<mnt6_swparams>(P, Q);
}

/* choice of pairing */
mnt6_G1_precomp mnt6_precompute_G1(const mnt6_G1& P)
{
    return short_weierstrass_precompute_G1<mnt6_swparams>(P);
}

mnt6_G2_precomp mnt6_precompute_G2(const mnt6_G2& Q)
{
    return short_weierstrass_precompute_G2<mnt6_swparams>(Q);
}

mnt6_Fq6 mnt6_miller_loop(const mnt6_G1_precomp &prec_P,
                          const mnt6_G2_precomp &prec_Q)
{
    return short_weierstrass_miller_loop<mnt6_swparams>(prec_P, prec_Q);
}

mnt6_Fq6 mnt6_double_miller_loop(const mnt6_G1_precomp &prec_P1,
                                 const mnt6_G2_precomp &prec_Q1,
                                 const mnt6_G1_precomp &prec_P2,
                                 const mnt6_G2_precomp &prec_Q2)
{
    return short_weierstrass_double_miller_loop<mnt6_swparams>(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

mnt6_Fq6 mnt6_pairing(const mnt6_G1& P,
                      const mnt6_G2 &Q)
{
    return short_weierstrass_pairing<mnt6_swparams>(P, Q);
}

mnt6_GT mnt6_reduced_pairing(const mnt6_G1 &P,
                             const mnt6_G2 &Q)
{
    return short_weierstrass_reduced_pairing<mnt6_swparams>(P, Q);
}

mnt6_GT mnt6_affine_reduced_pairing(const mnt6_G1 &P,
                                    const mnt6_G2 &Q)
{
    return short_weierstrass_affine_reduced_pairing<mnt6_swparams>(P, Q);
}


} // libff
