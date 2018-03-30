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

mnt6_Fq3 mnt6_swparams::mul_by_a(const mnt6_Fq3 &elt)
{
    return mnt6_Fq3(mnt6_twist_mul_by_a_c0 * elt.c1, mnt6_twist_mul_by_a_c1 * elt.c2, mnt6_twist_mul_by_a_c2 * elt.c0);
}

mnt6_Fq3 mnt6_swparams::mul_by_b(const mnt6_Fq3 &elt)
{
    return mnt6_Fq3(mnt6_twist_mul_by_b_c0 * elt.c0, mnt6_twist_mul_by_b_c1 * elt.c1, mnt6_twist_mul_by_b_c2 * elt.c2);
}

}
