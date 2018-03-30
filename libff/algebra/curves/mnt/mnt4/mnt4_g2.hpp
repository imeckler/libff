/** @file
 *****************************************************************************

 Declaration of interfaces for the MNT4 G2 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT4_G2_HPP_
#define MNT4_G2_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g2.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_init.hpp>

namespace libff {

typedef short_weierstrass_G2<mnt4_swparams> mnt4_G2;

} // libff

#endif // MNT4_G2_HPP_
