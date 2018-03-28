/** @file
 *****************************************************************************

 Declaration of interfaces for the MNT6 G1 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MNT6_G1_HPP_
#define MNT6_G1_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g1.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_init.hpp>

namespace libff {

typedef short_weierstrass_G1<mnt6_swparams> mnt6_G1;

} // libff

#endif // MNT6_G1_HPP_
