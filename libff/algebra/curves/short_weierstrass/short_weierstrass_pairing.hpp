/** @file
 *****************************************************************************

 Declaration of interfaces for pairing operations on short weierstrass curves

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SHORT_WEIERSTRASS_PAIRING_HPP_
#define SHORT_WEIERSTRASS_PAIRING_HPP_

#include <vector>

#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g1.hpp>
#include <libff/algebra/curves/short_weierstrass/short_weierstrass_g2.hpp>
#include <libff/common/profiling.hpp>

namespace libff {

/* final exponentiation */
template<typename SWParamsT>
using SWFqk = typename SWParamsT::Fqk;

template<typename SWParamsT>
using SWGT = typename SWParamsT::GT;

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_final_exponentiation_last_chunk(const SWFqk<SWParamsT> &elt,
                                                      const SWFqk<SWParamsT> &elt_inv);
template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_final_exponentiation(const SWFqk<SWParamsT> &elt);

/* affine ate miller loop */

template<typename SWParamsT>
struct short_weierstrass_affine_ate_G1_precomputation {
    SWFq<SWParamsT> PX;
    SWFq<SWParamsT> PY;
    SWtwist_field<SWParamsT> PY_twist_squared;
};

template<typename SWParamsT>
struct short_weierstrass_affine_ate_coeffs {
    // TODO: trim (not all of them are needed)
    SWtwist_field<SWParamsT> old_RX;
    SWtwist_field<SWParamsT> old_RY;
    SWtwist_field<SWParamsT> gamma;
    SWtwist_field<SWParamsT> gamma_twist;
    SWtwist_field<SWParamsT> gamma_X;
};

template<typename SWParamsT>
struct short_weierstrass_affine_ate_G2_precomputation {
    SWtwist_field<SWParamsT> QX;
    SWtwist_field<SWParamsT> QY;
    std::vector<short_weierstrass_affine_ate_coeffs<SWParamsT>> coeffs;
};

template<typename SWParamsT>
short_weierstrass_affine_ate_G1_precomputation<SWParamsT>
short_weierstrass_affine_ate_precompute_G1(const short_weierstrass_G1<SWParamsT>& P);
template<typename SWParamsT>
short_weierstrass_affine_ate_G2_precomputation<SWParamsT>
short_weierstrass_affine_ate_precompute_G2(const short_weierstrass_G2<SWParamsT>& Q);

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_affine_ate_miller_loop(
    const short_weierstrass_affine_ate_G1_precomputation<SWParamsT> &prec_P,
    const short_weierstrass_affine_ate_G2_precomputation<SWParamsT> &prec_Q);

/* ate pairing */

template<typename SWParamsT>
struct short_weierstrass_ate_G1_precomp {
    SWFq<SWParamsT> PX;
    SWFq<SWParamsT> PY;
    SWtwist_field<SWParamsT> PX_twist;
    SWtwist_field<SWParamsT> PY_twist;

    bool operator==(const short_weierstrass_ate_G1_precomp<SWParamsT> &other) const;
    template<typename SWParamsTT>
    friend std::ostream& operator<<(std::ostream &out, const short_weierstrass_ate_G1_precomp<SWParamsTT> &prec_P);
    template<typename SWParamsTT>
    friend std::istream& operator>>(std::istream &in, short_weierstrass_ate_G1_precomp<SWParamsTT> &prec_P);
};

template<typename SWParamsT>
struct short_weierstrass_ate_dbl_coeffs {
    SWtwist_field<SWParamsT> c_H;
    SWtwist_field<SWParamsT> c_4C;
    SWtwist_field<SWParamsT> c_J;
    SWtwist_field<SWParamsT> c_L;

    bool operator==(const short_weierstrass_ate_dbl_coeffs<SWParamsT> &other) const;
    template<typename SWParamsTT>
    friend std::ostream& operator<<(std::ostream &out, const short_weierstrass_ate_dbl_coeffs<SWParamsTT> &dc);
    template<typename SWParamsTT>
    friend std::istream& operator>>(std::istream &in, short_weierstrass_ate_dbl_coeffs<SWParamsTT> &dc);
};

template<typename SWParamsT>
struct short_weierstrass_ate_add_coeffs {
    SWtwist_field<SWParamsT> c_L1;
    SWtwist_field<SWParamsT> c_RZ;

    bool operator==(const short_weierstrass_ate_add_coeffs<SWParamsT> &other) const;
    template<typename SWParamsTT>
    friend std::ostream& operator<<(std::ostream &out, const short_weierstrass_ate_add_coeffs<SWParamsTT> &dc);
    template<typename SWParamsTT>
    friend std::istream& operator>>(std::istream &in, short_weierstrass_ate_add_coeffs<SWParamsTT> &dc);
};

template<typename SWParamsT>
struct short_weierstrass_ate_G2_precomp {
    SWtwist_field<SWParamsT> QX;
    SWtwist_field<SWParamsT> QY;
    SWtwist_field<SWParamsT> QY2;
    SWtwist_field<SWParamsT> QX_over_twist;
    SWtwist_field<SWParamsT> QY_over_twist;
    std::vector<short_weierstrass_ate_dbl_coeffs<SWParamsT>> dbl_coeffs;
    std::vector<short_weierstrass_ate_add_coeffs<SWParamsT>> add_coeffs;

    bool operator==(const short_weierstrass_ate_G2_precomp<SWParamsT> &other) const;
    template<typename SWParamsTT>
    friend std::ostream& operator<<(std::ostream &out, const short_weierstrass_ate_G2_precomp<SWParamsTT> &prec_Q);
    template<typename SWParamsTT>
    friend std::istream& operator>>(std::istream &in, short_weierstrass_ate_G2_precomp<SWParamsTT> &prec_Q);
};

template<typename SWParamsT>
short_weierstrass_ate_G1_precomp<SWParamsT>
short_weierstrass_ate_precompute_G1(const short_weierstrass_G1<SWParamsT>& P);
template<typename SWParamsT>
short_weierstrass_ate_G2_precomp<SWParamsT>
short_weierstrass_ate_precompute_G2(const short_weierstrass_G2<SWParamsT>& Q);

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_ate_miller_loop(
    const short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P,
    const short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q);
template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_ate_double_miller_loop(const short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P1,
                                                          const short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q1,
                                                          const short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P2,
                                                          const short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q2);

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_ate_pairing(
    const short_weierstrass_G1<SWParamsT>& P,
    const short_weierstrass_G2<SWParamsT> &Q);
template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_ate_reduced_pairing(
    const short_weierstrass_G1<SWParamsT> &P,
    const short_weierstrass_G2<SWParamsT> &Q);

/* choice of pairing */

template<typename SWParamsT>
using short_weierstrass_G1_precomp = short_weierstrass_ate_G1_precomp<SWParamsT>;

template<typename SWParamsT>
using short_weierstrass_G2_precomp = short_weierstrass_ate_G2_precomp<SWParamsT>;

template<typename SWParamsT>
short_weierstrass_G1_precomp<SWParamsT> short_weierstrass_precompute_G1(const short_weierstrass_G1<SWParamsT>& P);

template<typename SWParamsT>
short_weierstrass_G2_precomp<SWParamsT> short_weierstrass_precompute_G2(const short_weierstrass_G2<SWParamsT>& Q);

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_miller_loop(
    const short_weierstrass_G1_precomp<SWParamsT> &prec_P,
    const short_weierstrass_G2_precomp<SWParamsT> &prec_Q);

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_double_miller_loop(
    const short_weierstrass_G1_precomp<SWParamsT> &prec_P1,
    const short_weierstrass_G2_precomp<SWParamsT> &prec_Q1,
    const short_weierstrass_G1_precomp<SWParamsT> &prec_P2,
    const short_weierstrass_G2_precomp<SWParamsT> &prec_Q2);

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_pairing(
    const short_weierstrass_G1<SWParamsT>& P,
    const short_weierstrass_G2<SWParamsT> &Q);

template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_reduced_pairing(
    const short_weierstrass_G1<SWParamsT> &P,
    const short_weierstrass_G2<SWParamsT> &Q);

template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_affine_reduced_pairing(
    const short_weierstrass_G1<SWParamsT> &P,
    const short_weierstrass_G2<SWParamsT> &Q);

// Begin implementation
template<typename SWParamsT>
bool short_weierstrass_ate_G1_precomp<SWParamsT>::operator==(const short_weierstrass_ate_G1_precomp<SWParamsT> &other) const
{
    return (this->PX == other.PX &&
            this->PY == other.PY &&
            this->PX_twist == other.PX_twist &&
            this->PY_twist == other.PY_twist);
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream &out, const short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P)
{
    out << prec_P.PX << OUTPUT_SEPARATOR << prec_P.PY << OUTPUT_SEPARATOR << prec_P.PX_twist << OUTPUT_SEPARATOR << prec_P.PY_twist;

    return out;
}

template<typename SWParamsT>
std::istream& operator>>(std::istream &in, short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P)
{
    in >> prec_P.PX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PX_twist;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_P.PY_twist;

    return in;
}

template<typename SWParamsT>
bool short_weierstrass_ate_dbl_coeffs<SWParamsT>::operator==(const short_weierstrass_ate_dbl_coeffs<SWParamsT> &other) const
{
    return (this->c_H == other.c_H &&
            this->c_4C == other.c_4C &&
            this->c_J == other.c_J &&
            this->c_L == other.c_L);
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream &out, const short_weierstrass_ate_dbl_coeffs<SWParamsT> &dc)
{
    out << dc.c_H << OUTPUT_SEPARATOR << dc.c_4C << OUTPUT_SEPARATOR << dc.c_J << OUTPUT_SEPARATOR << dc.c_L;
    return out;
}

template<typename SWParamsT>
std::istream& operator>>(std::istream &in, short_weierstrass_ate_dbl_coeffs<SWParamsT> &dc)
{
    in >> dc.c_H;
    consume_OUTPUT_SEPARATOR(in);
    in >> dc.c_4C;
    consume_OUTPUT_SEPARATOR(in);
    in >> dc.c_J;
    consume_OUTPUT_SEPARATOR(in);
    in >> dc.c_L;

    return in;
}

template<typename SWParamsT>
bool short_weierstrass_ate_add_coeffs<SWParamsT>::operator==(const short_weierstrass_ate_add_coeffs<SWParamsT> &other) const
{
    return (this->c_L1 == other.c_L1 &&
            this->c_RZ == other.c_RZ);
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream &out, const short_weierstrass_ate_add_coeffs<SWParamsT> &ac)
{
    out << ac.c_L1 << OUTPUT_SEPARATOR << ac.c_RZ;
    return out;
}

template<typename SWParamsT>
std::istream& operator>>(std::istream &in, short_weierstrass_ate_add_coeffs<SWParamsT> &ac)
{
    in >> ac.c_L1;
    consume_OUTPUT_SEPARATOR(in);
    in >> ac.c_RZ;
    return in;
}

template<typename SWParamsT>
bool short_weierstrass_ate_G2_precomp<SWParamsT>::operator==(const short_weierstrass_ate_G2_precomp<SWParamsT> &other) const
{
    return (this->QX == other.QX &&
            this->QY == other.QY &&
            this->QY2 == other.QY2 &&
            this->QX_over_twist == other.QX_over_twist &&
            this->QY_over_twist == other.QY_over_twist &&
            this->dbl_coeffs == other.dbl_coeffs &&
            this->add_coeffs == other.add_coeffs);
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream& out, const short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q)
{
    out << prec_Q.QX << OUTPUT_SEPARATOR
        << prec_Q.QY << OUTPUT_SEPARATOR
        << prec_Q.QY2  << OUTPUT_SEPARATOR
        << prec_Q.QX_over_twist << OUTPUT_SEPARATOR
        << prec_Q.QY_over_twist << "\n";
    out << prec_Q.dbl_coeffs.size() << "\n";
    for (const short_weierstrass_ate_dbl_coeffs<SWParamsT> &dc : prec_Q.dbl_coeffs)
    {
        out << dc << OUTPUT_NEWLINE;
    }
    out << prec_Q.add_coeffs.size() << "\n";
    for (const short_weierstrass_ate_add_coeffs<SWParamsT> &ac : prec_Q.add_coeffs)
    {
        out << ac << OUTPUT_NEWLINE;
    }

    return out;
}

template<typename SWParamsT>
std::istream& operator>>(std::istream& in, short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q)
{
    in >> prec_Q.QX;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY2;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QX_over_twist;
    consume_OUTPUT_SEPARATOR(in);
    in >> prec_Q.QY_over_twist;
    consume_newline(in);

    prec_Q.dbl_coeffs.clear();
    size_t dbl_s;
    in >> dbl_s;
    consume_newline(in);

    prec_Q.dbl_coeffs.reserve(dbl_s);

    for (size_t i = 0; i < dbl_s; ++i)
    {
        short_weierstrass_ate_dbl_coeffs<SWParamsT> dc;
        in >> dc;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.dbl_coeffs.emplace_back(dc);
    }

    prec_Q.add_coeffs.clear();
    size_t add_s;
    in >> add_s;
    consume_newline(in);

    prec_Q.add_coeffs.reserve(add_s);

    for (size_t i = 0; i < add_s; ++i)
    {
        short_weierstrass_ate_add_coeffs<SWParamsT> ac;
        in >> ac;
        consume_OUTPUT_NEWLINE(in);
        prec_Q.add_coeffs.emplace_back(ac);
    }

    return in;
}

/* final exponentiations */

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_final_exponentiation_last_chunk(const SWFqk<SWParamsT> &elt, const SWFqk<SWParamsT> &elt_inv)
{
    enter_block("Call to mnt4_final_exponentiation_last_chunk");
    const SWFqk<SWParamsT> elt_q = elt.Frobenius_map(1);
    SWFqk<SWParamsT> w1_part = elt_q.cyclotomic_exp(SWParamsT::final_exponent_last_chunk_w1);
    SWFqk<SWParamsT> w0_part;
    if (SWParamsT::final_exponent_last_chunk_is_w0_neg)
    {
    	w0_part = elt_inv.cyclotomic_exp(SWParamsT::final_exponent_last_chunk_abs_of_w0);
    } else {
    	w0_part = elt.cyclotomic_exp(SWParamsT::final_exponent_last_chunk_abs_of_w0);
    }
    SWFqk<SWParamsT> result = w1_part * w0_part;
    leave_block("Call to mnt4_final_exponentiation_last_chunk");

    return result;
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_final_exponentiation_first_chunk(
    const SWFqk<SWParamsT> &elt, const SWFqk<SWParamsT> &elt_inv)
{
    enter_block("Call to mnt4_final_exponentiation_first_chunk");
    const SWFqk<SWParamsT> result = SWParamsT::final_exponentiation_first_chunk(elt, elt_inv);
    leave_block("Call to mnt4_final_exponentiation_first_chunk");

    return result;
}

template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_final_exponentiation(const SWFqk<SWParamsT> &elt)
{
    enter_block("Call to mnt4_final_exponentiation");
    const SWFqk<SWParamsT> elt_inv = elt.inverse();
    const SWFqk<SWParamsT> elt_to_first_chunk = short_weierstrass_final_exponentiation_first_chunk<SWParamsT>(elt, elt_inv);
    const SWFqk<SWParamsT> elt_inv_to_first_chunk = short_weierstrass_final_exponentiation_first_chunk<SWParamsT>(elt_inv, elt);
    SWGT<SWParamsT> result = short_weierstrass_final_exponentiation_last_chunk<SWParamsT>(elt_to_first_chunk, elt_inv_to_first_chunk);
    leave_block("Call to mnt4_final_exponentiation");

    return result;
}

/* affine ate miller loop */

template<typename SWParamsT>
short_weierstrass_affine_ate_G1_precomputation<SWParamsT> short_weierstrass_affine_ate_precompute_G1(const short_weierstrass_G1<SWParamsT>& P)
{
    enter_block("Call to mnt4_affine_ate_precompute_G1");

    short_weierstrass_G1<SWParamsT> Pcopy = P;
    Pcopy.to_affine_coordinates();

    short_weierstrass_affine_ate_G1_precomputation<SWParamsT> result;
    result.PX = Pcopy.X();
    result.PY = Pcopy.Y();
    result.PY_twist_squared = Pcopy.Y() * SWParamsT::twist.squared();

    leave_block("Call to mnt4_affine_ate_precompute_G1");
    return result;
}

template<typename SWParamsT>
short_weierstrass_affine_ate_G2_precomputation<SWParamsT> short_weierstrass_affine_ate_precompute_G2(const short_weierstrass_G2<SWParamsT>& Q)
{
    enter_block("Call to mnt4_affine_ate_precompute_G2");

    short_weierstrass_G2<SWParamsT> Qcopy(Q);
    Qcopy.to_affine_coordinates();

    short_weierstrass_affine_ate_G2_precomputation<SWParamsT> result;
    result.QX = Qcopy.X();
    result.QY = Qcopy.Y();

    SWtwist_field<SWParamsT> RX = Qcopy.X();
    SWtwist_field<SWParamsT> RY = Qcopy.Y();

    const bigint<SWFr<SWParamsT>::num_limbs> &loop_count = SWParamsT::ate_loop_count;
    bool found_nonzero = false;

    std::vector<long> NAF = find_wnaf(1, loop_count);
    for (long i = NAF.size() - 1; i >= 0; --i)
    {
        if (!found_nonzero)
        {
            /* this skips the MSB itself */
            found_nonzero |= (NAF[i] != 0);
            continue;
        }

        short_weierstrass_affine_ate_coeffs<SWParamsT> c;
        c.old_RX = RX;
        c.old_RY = RY;
        SWtwist_field<SWParamsT> old_RX_2 = c.old_RX.squared();
        c.gamma = (old_RX_2 + old_RX_2 + old_RX_2 + SWParamsT::twist_coeff_a) * (c.old_RY + c.old_RY).inverse();
        c.gamma_twist = c.gamma * SWParamsT::twist;
        c.gamma_X = c.gamma * c.old_RX;
        result.coeffs.push_back(c);

        RX = c.gamma.squared() - (c.old_RX+c.old_RX);
        RY = c.gamma * (c.old_RX - RX) - c.old_RY;

        if (NAF[i] != 0)
        {
            short_weierstrass_affine_ate_coeffs<SWParamsT> c;
            c.old_RX = RX;
            c.old_RY = RY;
            if (NAF[i] > 0)
            {
                c.gamma = (c.old_RY - result.QY) * (c.old_RX - result.QX).inverse();
            }
            else
            {
                c.gamma = (c.old_RY + result.QY) * (c.old_RX - result.QX).inverse();
            }
            c.gamma_twist = c.gamma * SWParamsT::twist;
            c.gamma_X = c.gamma * result.QX;
            result.coeffs.push_back(c);

            RX = c.gamma.squared() - (c.old_RX+result.QX);
            RY = c.gamma * (c.old_RX - RX) - c.old_RY;
        }
    }

    /* TODO: maybe handle neg
       if (mnt4_ate_is_loop_count_neg)
       {
       mnt4_ate_add_coeffs ac;
       mnt4_affine_ate_dbl_coeffs c;
       c.old_RX = RX;
       c.old_RY = -RY;
       old_RX_2 = c.old_RY.squared();
       c.gamma = (old_RX_2 + old_RX_2 + old_RX_2 + mnt4_coeff_a) * (c.old_RY + c.old_RY).inverse();
       c.gamma_twist = c.gamma * mnt4_twist;
       c.gamma_X = c.gamma * c.old_RX;
       result.coeffs.push_back(c);
       }
    */

    leave_block("Call to mnt4_affine_ate_precompute_G2");
    return result;
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_affine_ate_miller_loop(const short_weierstrass_affine_ate_G1_precomputation<SWParamsT> &prec_P,
                                             const short_weierstrass_affine_ate_G2_precomputation<SWParamsT> &prec_Q)
{
    enter_block("Call to mnt4_affine_ate_miller_loop");

    SWFqk<SWParamsT> f = SWFqk<SWParamsT>::one();

    bool found_nonzero = false;
    size_t idx = 0;
    const bigint<SWFr<SWParamsT>::num_limbs> &loop_count = SWParamsT::ate_loop_count;

    std::vector<long> NAF = find_wnaf(1, loop_count);
    for (long i = NAF.size() - 1; i >= 0; --i)
    {
        if (!found_nonzero)
        {
            /* this skips the MSB itself */
            found_nonzero |= (NAF[i] != 0);
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           mnt4_param_p (skipping leading zeros) in MSB to LSB
           order */
        short_weierstrass_affine_ate_coeffs<SWParamsT> c = prec_Q.coeffs[idx++];

        SWFqk<SWParamsT> g_RR_at_P = SWFqk<SWParamsT>(prec_P.PY_twist_squared,
                                      - prec_P.PX * c.gamma_twist + c.gamma_X - c.old_RY);
        f = SWParamsT::special_mul(f.squared(), g_RR_at_P);

        if (NAF[i] != 0)
        {
            short_weierstrass_affine_ate_coeffs<SWParamsT> c = prec_Q.coeffs[idx++];
            SWFqk<SWParamsT> g_RQ_at_P;
            if (NAF[i] > 0)
            {
                g_RQ_at_P = SWFqk<SWParamsT>(prec_P.PY_twist_squared,
                                     - prec_P.PX * c.gamma_twist + c.gamma_X - prec_Q.QY);
            }
            else
            {
                g_RQ_at_P = SWFqk<SWParamsT>(prec_P.PY_twist_squared,
                                     - prec_P.PX * c.gamma_twist + c.gamma_X + prec_Q.QY);
            }
            f = SWParamsT::special_mul(f, g_RQ_at_P);
        }
    }

    /* TODO: maybe handle neg
       if (mnt4_ate_is_loop_count_neg)
       {
       // TODO:
       mnt4_affine_ate_coeffs ac = prec_Q.coeffs[idx++];
       SWFqk<SWParamsT> g_RnegR_at_P = SWFqk<SWParamsT>(prec_P.PY_twist_squared,
       - prec_P.PX * c.gamma_twist + c.gamma_X - c.old_RY);
       f = (f * g_RnegR_at_P).inverse();
       }
    */

    leave_block("Call to mnt4_affine_ate_miller_loop");

    return f;
}

/* ate pairing */

template<typename SWParamsT>
struct extended_short_weierstrass_G2_projective {
    SWtwist_field<SWParamsT> X;
    SWtwist_field<SWParamsT> Y;
    SWtwist_field<SWParamsT> Z;
    SWtwist_field<SWParamsT> T;

    void print() const
    {
        printf("extended short_weierstrass_G2<SWParamsT> projective X/Y/Z/T:\n");
        X.print();
        Y.print();
        Z.print();
        T.print();
    }

    void test_invariant() const
    {
        assert(T == Z.squared());
    }
};

template<typename SWParamsT>
void doubling_step_for_flipped_miller_loop(
    extended_short_weierstrass_G2_projective<SWParamsT> &current,
    short_weierstrass_ate_dbl_coeffs<SWParamsT> &dc)
{
    const SWtwist_field<SWParamsT> X = current.X, Y = current.Y, Z = current.Z, T = current.T;

    const SWtwist_field<SWParamsT> A = T.squared(); // A = T1^2
    const SWtwist_field<SWParamsT> B = X.squared(); // B = X1^2
    const SWtwist_field<SWParamsT> C = Y.squared(); // C = Y1^2
    const SWtwist_field<SWParamsT> D = C.squared(); // D = C^2
    const SWtwist_field<SWParamsT> E = (X+C).squared() - B - D; // E = (X1+C)^2-B-D
    const SWtwist_field<SWParamsT> F = (B+B+B) + SWParamsT::twist_coeff_a * A; // F = 3*B +  a  *A
    const SWtwist_field<SWParamsT> G = F.squared(); // G = F^2

    current.X = -(E+E+E+E) + G; // X3 = -4*E+G
    current.Y = -SWFq<SWParamsT>("8")*D + F*(E+E-current.X); // Y3 = -8*D+F*(2*E-X3)
    current.Z = (Y+Z).squared() - C - Z.squared(); // Z3 = (Y1+Z1)^2-C-Z1^2
    current.T = current.Z.squared(); // T3 = Z3^2

    dc.c_H = (current.Z + T).squared() - current.T - A; // H = (Z3+T1)^2-T3-A
    dc.c_4C = C+C+C+C; // fourC = 4*C
    dc.c_J = (F+T).squared() - G - A; // J = (F+T1)^2-G-A
    dc.c_L = (F+X).squared() - G - B; // L = (F+X1)^2-G-B

#ifdef DEBUG
    current.test_invariant();
#endif
}

template<typename SWParamsT>
void mixed_addition_step_for_flipped_miller_loop(
    const SWtwist_field<SWParamsT> base_X, const SWtwist_field<SWParamsT> base_Y, const SWtwist_field<SWParamsT> base_Y_squared,
    extended_short_weierstrass_G2_projective<SWParamsT> &current,
    short_weierstrass_ate_add_coeffs<SWParamsT> &ac)
{
    const SWtwist_field<SWParamsT> X1 = current.X, Y1 = current.Y, Z1 = current.Z, T1 = current.T;
    const SWtwist_field<SWParamsT> &x2 = base_X,    &y2 =  base_Y, &y2_squared = base_Y_squared;

    const SWtwist_field<SWParamsT> B = x2 * T1; // B = x2 * T1
    const SWtwist_field<SWParamsT> D = ((y2 + Z1).squared() - y2_squared - T1) * T1; // D = ((y2 + Z1)^2 - y2squared - T1) * T1
    const SWtwist_field<SWParamsT> H = B - X1; // H = B - X1
    const SWtwist_field<SWParamsT> I = H.squared(); // I = H^2
    const SWtwist_field<SWParamsT> E = I + I + I + I; // E = 4*I
    const SWtwist_field<SWParamsT> J = H * E; // J = H * E
    const SWtwist_field<SWParamsT> V = X1 * E; // V = X1 * E
    const SWtwist_field<SWParamsT> L1 = D - (Y1 + Y1); // L1 = D - 2 * Y1

    current.X = L1.squared() - J - (V+V); // X3 = L1^2 - J - 2*V
    current.Y = L1 * (V-current.X) - (Y1+Y1) * J; // Y3 = L1 * (V-X3) - 2*Y1 * J
    current.Z = (Z1+H).squared() - T1 - I; // Z3 = (Z1 + H)^2 - T1 - I
    current.T = current.Z.squared(); // T3 = Z3^2

    ac.c_L1 = L1;
    ac.c_RZ = current.Z;
#ifdef DEBUG
    current.test_invariant();
#endif
}

template<typename SWParamsT>
short_weierstrass_ate_G1_precomp<SWParamsT> short_weierstrass_ate_precompute_G1(const short_weierstrass_G1<SWParamsT>& P)
{
    enter_block("Call to mnt4_ate_precompute_G1");

    short_weierstrass_G1<SWParamsT> Pcopy = P;
    Pcopy.to_affine_coordinates();

    short_weierstrass_ate_G1_precomp<SWParamsT> result;
    result.PX = Pcopy.X();
    result.PY = Pcopy.Y();
    result.PX_twist = Pcopy.X() * SWParamsT::twist;
    result.PY_twist = Pcopy.Y() * SWParamsT::twist;

    leave_block("Call to mnt4_ate_precompute_G1");
    return result;
}

template<typename SWParamsT>
short_weierstrass_ate_G2_precomp<SWParamsT> short_weierstrass_ate_precompute_G2(const short_weierstrass_G2<SWParamsT>& Q)
{
    enter_block("Call to mnt4_ate_precompute_G2");

    short_weierstrass_G2<SWParamsT> Qcopy(Q);
    Qcopy.to_affine_coordinates();

    short_weierstrass_ate_G2_precomp<SWParamsT> result;
    result.QX = Qcopy.X();
    result.QY = Qcopy.Y();
    result.QY2 = Qcopy.Y().squared();
    result.QX_over_twist = Qcopy.X() * SWParamsT::twist.inverse();
    result.QY_over_twist = Qcopy.Y() * SWParamsT::twist.inverse();

    extended_short_weierstrass_G2_projective<SWParamsT> R;
    R.X = Qcopy.X();
    R.Y = Qcopy.Y();
    R.Z = SWtwist_field<SWParamsT>::one();
    R.T = SWtwist_field<SWParamsT>::one();

    const bigint<SWFr<SWParamsT>::num_limbs> &loop_count = SWParamsT::ate_loop_count;
    bool found_one = false;

    for (long i = loop_count.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);
        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        short_weierstrass_ate_dbl_coeffs<SWParamsT> dc;
        doubling_step_for_flipped_miller_loop(R, dc);
        result.dbl_coeffs.push_back(dc);
        if (bit)
        {
            short_weierstrass_ate_add_coeffs<SWParamsT> ac;
            mixed_addition_step_for_flipped_miller_loop(result.QX, result.QY, result.QY2, R, ac);
            result.add_coeffs.push_back(ac);
        }
    }

    if (SWParamsT::ate_is_loop_count_neg)
    {
    	SWtwist_field<SWParamsT> RZ_inv = R.Z.inverse();
    	SWtwist_field<SWParamsT> RZ2_inv = RZ_inv.squared();
    	SWtwist_field<SWParamsT> RZ3_inv = RZ2_inv * RZ_inv;
    	SWtwist_field<SWParamsT> minus_R_affine_X = R.X * RZ2_inv;
    	SWtwist_field<SWParamsT> minus_R_affine_Y = - R.Y * RZ3_inv;
    	SWtwist_field<SWParamsT> minus_R_affine_Y2 = minus_R_affine_Y.squared();
    	short_weierstrass_ate_add_coeffs<SWParamsT> ac;
        mixed_addition_step_for_flipped_miller_loop(minus_R_affine_X, minus_R_affine_Y, minus_R_affine_Y2, R, ac);
        result.add_coeffs.push_back(ac);
    }

    leave_block("Call to mnt4_ate_precompute_G2");
    return result;
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_ate_miller_loop(
    const short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P,
    const short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q)
{
    enter_block("Call to mnt4_ate_miller_loop");

    SWtwist_field<SWParamsT> L1_coeff = SWParamsT::embed_as_first_coordinate(prec_P.PX) - prec_Q.QX_over_twist;

    SWFqk<SWParamsT> f = SWFqk<SWParamsT>::one();

    bool found_one = false;
    size_t dbl_idx = 0;
    size_t add_idx = 0;

    const bigint<SWFr<SWParamsT>::num_limbs> &loop_count = SWParamsT::ate_loop_count;
    for (long i = loop_count.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);

        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           mnt4_param_p (skipping leading zeros) in MSB to LSB
           order */
        short_weierstrass_ate_dbl_coeffs<SWParamsT> dc = prec_Q.dbl_coeffs[dbl_idx++];

        SWFqk<SWParamsT> g_RR_at_P = SWFqk<SWParamsT>(- dc.c_4C - dc.c_J * prec_P.PX_twist + dc.c_L,
                                      dc.c_H * prec_P.PY_twist);
        f = f.squared() * g_RR_at_P;
        if (bit)
        {
            short_weierstrass_ate_add_coeffs<SWParamsT> ac = prec_Q.add_coeffs[add_idx++];

            SWFqk<SWParamsT> g_RQ_at_P = SWFqk<SWParamsT>(ac.c_RZ * prec_P.PY_twist,
                                          -(prec_Q.QY_over_twist * ac.c_RZ + L1_coeff * ac.c_L1));
            f = f * g_RQ_at_P;
        }
    }

    if (SWParamsT::ate_is_loop_count_neg)
    {
    	short_weierstrass_ate_add_coeffs<SWParamsT> ac = prec_Q.add_coeffs[add_idx++];
    	SWFqk<SWParamsT> g_RnegR_at_P = SWFqk<SWParamsT>(ac.c_RZ * prec_P.PY_twist,
                                         -(prec_Q.QY_over_twist * ac.c_RZ + L1_coeff * ac.c_L1));
    	f = (f * g_RnegR_at_P).inverse();
    }

    leave_block("Call to mnt4_ate_miller_loop");

    return f;
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_ate_double_miller_loop(
    const short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P1,
    const short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q1,
    const short_weierstrass_ate_G1_precomp<SWParamsT> &prec_P2,
    const short_weierstrass_ate_G2_precomp<SWParamsT> &prec_Q2)
{
    enter_block("Call to mnt4_ate_double_miller_loop");

    SWtwist_field<SWParamsT> L1_coeff1 = SWParamsT::embed_as_first_coordinate(prec_P1.PX) - prec_Q1.QX_over_twist;
    SWtwist_field<SWParamsT> L1_coeff2 = SWParamsT::embed_as_first_coordinate(prec_P2.PX) - prec_Q2.QX_over_twist;

    SWFqk<SWParamsT> f = SWFqk<SWParamsT>::one();

    bool found_one = false;
    size_t dbl_idx = 0;
    size_t add_idx = 0;

    const bigint<SWFr<SWParamsT>::num_limbs> &loop_count = SWParamsT::ate_loop_count;
    for (long i = loop_count.max_bits() - 1; i >= 0; --i)
    {
        const bool bit = loop_count.test_bit(i);

        if (!found_one)
        {
            /* this skips the MSB itself */
            found_one |= bit;
            continue;
        }

        /* code below gets executed for all bits (EXCEPT the MSB itself) of
           mnt4_param_p (skipping leading zeros) in MSB to LSB
           order */
        short_weierstrass_ate_dbl_coeffs<SWParamsT> dc1 = prec_Q1.dbl_coeffs[dbl_idx];
        short_weierstrass_ate_dbl_coeffs<SWParamsT> dc2 = prec_Q2.dbl_coeffs[dbl_idx];
        ++dbl_idx;

        SWFqk<SWParamsT> g_RR_at_P1 = SWFqk<SWParamsT>(- dc1.c_4C - dc1.c_J * prec_P1.PX_twist + dc1.c_L,
                                       dc1.c_H * prec_P1.PY_twist);

        SWFqk<SWParamsT> g_RR_at_P2 = SWFqk<SWParamsT>(- dc2.c_4C - dc2.c_J * prec_P2.PX_twist + dc2.c_L,
                                       dc2.c_H * prec_P2.PY_twist);

        f = f.squared() * g_RR_at_P1 * g_RR_at_P2;

        if (bit)
        {
            short_weierstrass_ate_add_coeffs<SWParamsT> ac1 = prec_Q1.add_coeffs[add_idx];
            short_weierstrass_ate_add_coeffs<SWParamsT> ac2 = prec_Q2.add_coeffs[add_idx];
            ++add_idx;

            SWFqk<SWParamsT> g_RQ_at_P1 = SWFqk<SWParamsT>(ac1.c_RZ * prec_P1.PY_twist,
                                           -(prec_Q1.QY_over_twist * ac1.c_RZ + L1_coeff1 * ac1.c_L1));
            SWFqk<SWParamsT> g_RQ_at_P2 = SWFqk<SWParamsT>(ac2.c_RZ * prec_P2.PY_twist,
                                           -(prec_Q2.QY_over_twist * ac2.c_RZ + L1_coeff2 * ac2.c_L1));

            f = f * g_RQ_at_P1 * g_RQ_at_P2;
        }
    }

    if (SWParamsT::ate_is_loop_count_neg)
    {
    	short_weierstrass_ate_add_coeffs<SWParamsT> ac1 = prec_Q1.add_coeffs[add_idx];
        short_weierstrass_ate_add_coeffs<SWParamsT> ac2 = prec_Q2.add_coeffs[add_idx];
    	++add_idx;
    	SWFqk<SWParamsT> g_RnegR_at_P1 = SWFqk<SWParamsT>(ac1.c_RZ * prec_P1.PY_twist,
                                          -(prec_Q1.QY_over_twist * ac1.c_RZ + L1_coeff1 * ac1.c_L1));
    	SWFqk<SWParamsT> g_RnegR_at_P2 = SWFqk<SWParamsT>(ac2.c_RZ * prec_P2.PY_twist,
                                          -(prec_Q2.QY_over_twist * ac2.c_RZ + L1_coeff2 * ac2.c_L1));

    	f = (f * g_RnegR_at_P1 * g_RnegR_at_P2).inverse();
    }

    leave_block("Call to mnt4_ate_double_miller_loop");

    return f;
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_ate_pairing(const short_weierstrass_G1<SWParamsT>& P, const short_weierstrass_G2<SWParamsT> &Q)
{
    enter_block("Call to mnt4_ate_pairing");
    short_weierstrass_ate_G1_precomp<SWParamsT> prec_P = short_weierstrass_ate_precompute_G1<SWParamsT>(P);
    short_weierstrass_ate_G2_precomp<SWParamsT> prec_Q = short_weierstrass_ate_precompute_G2<SWParamsT>(Q);
    SWFqk<SWParamsT> result = short_weierstrass_ate_miller_loop<SWParamsT>(prec_P, prec_Q);
    leave_block("Call to mnt4_ate_pairing");
    return result;
}

template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_ate_reduced_pairing(const short_weierstrass_G1<SWParamsT> &P, const short_weierstrass_G2<SWParamsT> &Q)
{
    enter_block("Call to mnt4_ate_reduced_pairing");
    const SWFqk<SWParamsT> f = short_weierstrass_ate_pairing<SWParamsT>(P, Q);
    const SWGT<SWParamsT> result = short_weierstrass_final_exponentiation<SWParamsT>(f);
    leave_block("Call to mnt4_ate_reduced_pairing");
    return result;
}

template<typename SWParamsT>
short_weierstrass_G1_precomp<SWParamsT> short_weierstrass_precompute_G1(const short_weierstrass_G1<SWParamsT>& P)
{
    return short_weierstrass_ate_precompute_G1<SWParamsT>(P);
}

template<typename SWParamsT>
short_weierstrass_G2_precomp<SWParamsT> short_weierstrass_precompute_G2(const short_weierstrass_G2<SWParamsT>& Q)
{
    return short_weierstrass_ate_precompute_G2<SWParamsT>(Q);
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_miller_loop(
    const short_weierstrass_G1_precomp<SWParamsT> &prec_P,
    const short_weierstrass_G2_precomp<SWParamsT> &prec_Q)
{
    return short_weierstrass_ate_miller_loop<SWParamsT>(prec_P, prec_Q);
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_double_miller_loop(
    const short_weierstrass_G1_precomp<SWParamsT> &prec_P1,
    const short_weierstrass_G2_precomp<SWParamsT> &prec_Q1,
    const short_weierstrass_G1_precomp<SWParamsT> &prec_P2,
    const short_weierstrass_G2_precomp<SWParamsT> &prec_Q2)
{
    return short_weierstrass_ate_double_miller_loop<SWParamsT>(prec_P1, prec_Q1, prec_P2, prec_Q2);
}

template<typename SWParamsT>
SWFqk<SWParamsT> short_weierstrass_pairing(
    const short_weierstrass_G1<SWParamsT>& P,
    const short_weierstrass_G2<SWParamsT> &Q)
{
    return short_weierstrass_ate_pairing<SWParamsT>(P, Q);
}

template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_reduced_pairing(
    const short_weierstrass_G1<SWParamsT> &P,
    const short_weierstrass_G2<SWParamsT> &Q)
{
    return short_weierstrass_ate_reduced_pairing<SWParamsT>(P, Q);
}

template<typename SWParamsT>
SWGT<SWParamsT> short_weierstrass_affine_reduced_pairing(
    const short_weierstrass_G1<SWParamsT> &P,
    const short_weierstrass_G2<SWParamsT> &Q)
{
    const short_weierstrass_affine_ate_G1_precomputation<SWParamsT> prec_P = short_weierstrass_affine_ate_precompute_G1<SWParamsT>(P);
    const short_weierstrass_affine_ate_G2_precomputation<SWParamsT> prec_Q = short_weierstrass_affine_ate_precompute_G2<SWParamsT>(Q);
    const SWFqk<SWParamsT> f = short_weierstrass_affine_ate_miller_loop<SWParamsT>(prec_P, prec_Q);
    const SWGT<SWParamsT> result = short_weierstrass_final_exponentiation<SWParamsT>(f);
    return result;
}

} // libff

#endif // SHORT_WEIERSTRASS_PAIRING_HPP_
