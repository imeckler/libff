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
mnt4_Fq mnt4_swparams::zeroX;
mnt4_Fq mnt4_swparams::zeroY;
mnt4_Fq mnt4_swparams::zeroZ;
mnt4_Fq mnt4_swparams::oneX;
mnt4_Fq mnt4_swparams::oneY;
mnt4_Fq mnt4_swparams::oneZ;

void init_mnt4_swparams()
{
    /* choice of short Weierstrass curve */
    mnt4_swparams::coeff_a = mnt4_Fq("2");
    mnt4_swparams::coeff_b = mnt4_Fq("423894536526684178289416011533888240029318103673896002803341544124054745019340795360841685");
    mnt4_swparams::zeroX = mnt4_Fq::zero();
    mnt4_swparams::zeroY = mnt4_Fq::one();
    mnt4_swparams::zeroZ = mnt4_Fq::zero();
    mnt4_swparams::oneX = mnt4_Fq("60760244141852568949126569781626075788424196370144486719385562369396875346601926534016838");
    mnt4_swparams::oneY = mnt4_Fq("363732850702582978263902770815145784459747722357071843971107674179038674942891694705904306");
    mnt4_swparams::oneZ = mnt4_Fq::one();
}

}
