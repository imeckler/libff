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

mnt4_Fq2 mnt4_swparams::twist_coeff_a;
mnt4_Fq2 mnt4_swparams::twist_coeff_b;

mnt4_Fq2 mnt4_swparams::G2_zero_X;
mnt4_Fq2 mnt4_swparams::G2_zero_Y;
mnt4_Fq2 mnt4_swparams::G2_zero_Z;
mnt4_Fq2 mnt4_swparams::G2_one_X;
mnt4_Fq2 mnt4_swparams::G2_one_Y;
mnt4_Fq2 mnt4_swparams::G2_one_Z;

mnt4_Fq2 mnt4_swparams::twist_mul_by_q_X;
mnt4_Fq2 mnt4_swparams::twist_mul_by_q_Y;

mnt4_Fq2 mnt4_swparams::mul_by_a(const mnt4_Fq2 &elt)
{
    return mnt4_Fq2(mnt4_twist_mul_by_a_c0 * elt.c0, mnt4_twist_mul_by_a_c1 * elt.c1);
}

mnt4_Fq2 mnt4_swparams::mul_by_b(const mnt4_Fq2 &elt)
{
    return mnt4_Fq2(mnt4_twist_mul_by_b_c0 * elt.c1, mnt4_twist_mul_by_b_c1 * elt.c0);
}

}
