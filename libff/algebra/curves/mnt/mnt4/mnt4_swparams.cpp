/** @file
 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>

namespace libff {

mnt4_Fq mnt4_swparams::coeff_a;
mnt4_Fq mnt4_swparams::coeff_b;
mnt4_Fq mnt4_swparams::G1_zero_X;
mnt4_Fq mnt4_swparams::G1_zero_Y;
mnt4_Fq mnt4_swparams::G1_zero_Z;
mnt4_Fq mnt4_swparams::G1_one_X;
mnt4_Fq mnt4_swparams::G1_one_Y;
mnt4_Fq mnt4_swparams::G1_one_Z;

mnt4_Fq2 mnt4_swparams::twist;
mnt4_Fq2 mnt4_swparams::twist_coeff_a;
mnt4_Fq2 mnt4_swparams::twist_coeff_b;

mnt4_Fq2 mnt4_swparams::G2_zero_X;
mnt4_Fq2 mnt4_swparams::G2_zero_Y;
mnt4_Fq2 mnt4_swparams::G2_zero_Z;
mnt4_Fq2 mnt4_swparams::G2_one_X;
mnt4_Fq2 mnt4_swparams::G2_one_Y;
mnt4_Fq2 mnt4_swparams::G2_one_Z;

mnt4_Fq mnt4_swparams::twist_mul_by_q_X;
mnt4_Fq mnt4_swparams::twist_mul_by_q_Y;

bigint<mnt4_q_limbs> mnt4_swparams::ate_loop_count;
bool mnt4_swparams::ate_is_loop_count_neg;
bigint<4*mnt4_q_limbs> mnt4_swparams::final_exponent;
bigint<mnt4_q_limbs> mnt4_swparams::final_exponent_last_chunk_abs_of_w0;
bool mnt4_swparams::final_exponent_last_chunk_is_w0_neg;
bigint<mnt4_q_limbs> mnt4_swparams::final_exponent_last_chunk_w1;

mnt4_Fq2 mnt4_swparams::mul_by_a(const mnt4_Fq2 &elt)
{
    return mnt4_Fq2(mnt4_twist_mul_by_a_c0 * elt.c0, mnt4_twist_mul_by_a_c1 * elt.c1);
}

mnt4_Fq2 mnt4_swparams::mul_by_b(const mnt4_Fq2 &elt)
{
    return mnt4_Fq2(mnt4_twist_mul_by_b_c0 * elt.c1, mnt4_twist_mul_by_b_c1 * elt.c0);
}


mnt4_Fq4 mnt4_swparams::final_exponentiation_first_chunk(const mnt4_Fq4 &elt, const mnt4_Fq4 &elt_inv)
{
    /* (q^2-1) */

    /* elt_q2 = elt^(q^2) */
    const mnt4_Fq4 elt_q2 = elt.Frobenius_map(2);
    /* elt_q3_over_elt = elt^(q^2-1) */
    const mnt4_Fq4 elt_q2_over_elt = elt_q2 * elt_inv;
    return elt_q2_over_elt;
}

mnt4_Fq4 mnt4_swparams::special_mul(const mnt4_Fq4 &x, const mnt4_Fq4 &y)
{
    return x.mul_by_023(y);
}

mnt4_Fq2 mnt4_swparams::embed_as_first_coordinate(const mnt4_Fq &x)
{
    return mnt4_Fq2(x, mnt4_Fq::zero());
}

}
