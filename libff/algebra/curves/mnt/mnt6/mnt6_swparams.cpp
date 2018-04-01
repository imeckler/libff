/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/mnt/mnt6/mnt6_init.hpp>

namespace libff {

mnt6_Fq mnt6_swparams::coeff_a;
mnt6_Fq mnt6_swparams::coeff_b;
mnt6_Fq mnt6_swparams::G1_zero_X;
mnt6_Fq mnt6_swparams::G1_zero_Y;
mnt6_Fq mnt6_swparams::G1_zero_Z;
mnt6_Fq mnt6_swparams::G1_one_X;
mnt6_Fq mnt6_swparams::G1_one_Y;
mnt6_Fq mnt6_swparams::G1_one_Z;

mnt6_Fq3 mnt6_swparams::twist;
mnt6_Fq3 mnt6_swparams::twist_coeff_a;
mnt6_Fq3 mnt6_swparams::twist_coeff_b;

mnt6_Fq3 mnt6_swparams::G2_zero_X;
mnt6_Fq3 mnt6_swparams::G2_zero_Y;
mnt6_Fq3 mnt6_swparams::G2_zero_Z;
mnt6_Fq3 mnt6_swparams::G2_one_X;
mnt6_Fq3 mnt6_swparams::G2_one_Y;
mnt6_Fq3 mnt6_swparams::G2_one_Z;

mnt6_Fq mnt6_swparams::twist_mul_by_q_X;
mnt6_Fq mnt6_swparams::twist_mul_by_q_Y;

bigint<mnt6_q_limbs> mnt6_swparams::ate_loop_count;
bool mnt6_swparams::ate_is_loop_count_neg;
bigint<6*mnt6_q_limbs> mnt6_swparams::final_exponent;
bigint<mnt6_q_limbs> mnt6_swparams::final_exponent_last_chunk_abs_of_w0;
bool mnt6_swparams::final_exponent_last_chunk_is_w0_neg;
bigint<mnt6_q_limbs> mnt6_swparams::final_exponent_last_chunk_w1;

mnt6_Fq3 mnt6_swparams::mul_by_a(const mnt6_Fq3 &elt)
{
    return mnt6_Fq3(mnt6_twist_mul_by_a_c0 * elt.c1, mnt6_twist_mul_by_a_c1 * elt.c2, mnt6_twist_mul_by_a_c2 * elt.c0);
}

mnt6_Fq3 mnt6_swparams::mul_by_b(const mnt6_Fq3 &elt)
{
    return mnt6_Fq3(mnt6_twist_mul_by_b_c0 * elt.c0, mnt6_twist_mul_by_b_c1 * elt.c1, mnt6_twist_mul_by_b_c2 * elt.c2);
}

mnt6_Fq6 mnt6_swparams::final_exponentiation_first_chunk(const mnt6_Fq6 &elt, const mnt6_Fq6 &elt_inv)
{
    /* (q^3-1)*(q+1) */

    /* elt_q3 = elt^(q^3) */
    const mnt6_Fq6 elt_q3 = elt.Frobenius_map(3);
    /* elt_q3_over_elt = elt^(q^3-1) */
    const mnt6_Fq6 elt_q3_over_elt = elt_q3 * elt_inv;
    /* alpha = elt^((q^3-1) * q) */
    const mnt6_Fq6 alpha = elt_q3_over_elt.Frobenius_map(1);
    /* beta = elt^((q^3-1)*(q+1) */
    const mnt6_Fq6 beta = alpha * elt_q3_over_elt;
    return beta;
}

mnt6_Fq6 mnt6_swparams::special_mul(const mnt6_Fq6 &x, const mnt6_Fq6 &y)
{
    return x.mul_by_2345(y);
}

mnt6_Fq3 mnt6_swparams::embed_as_first_coordinate(const mnt6_Fq &x)
{
    return mnt6_Fq3(x, mnt6_Fq::zero(), mnt6_Fq::zero());
}

}
