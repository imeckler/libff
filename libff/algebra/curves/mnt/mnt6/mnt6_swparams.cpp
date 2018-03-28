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
mnt6_Fq mnt6_swparams::zeroX;
mnt6_Fq mnt6_swparams::zeroY;
mnt6_Fq mnt6_swparams::zeroZ;
mnt6_Fq mnt6_swparams::oneX;
mnt6_Fq mnt6_swparams::oneY;
mnt6_Fq mnt6_swparams::oneZ;

void init_mnt6_swparams()
{
    mnt6_swparams::coeff_a = mnt6_Fq("11");
    mnt6_swparams::coeff_b = mnt6_Fq("106700080510851735677967319632585352256454251201367587890185989362936000262606668469523074");
    mnt6_swparams::zeroX = mnt6_Fq::zero();
    mnt6_swparams::zeroY = mnt6_Fq::one();
    mnt6_swparams::zeroZ = mnt6_Fq::zero();
    mnt6_swparams::oneX = mnt6_Fq("336685752883082228109289846353937104185698209371404178342968838739115829740084426881123453");
    mnt6_swparams::oneY = mnt6_Fq("402596290139780989709332707716568920777622032073762749862342374583908837063963736098549800");
    mnt6_swparams::oneZ = mnt6_Fq::one();
}

}
