/** @file
 *****************************************************************************

 Declaration of interfaces for a twisted short weierstrass group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SHORT_WEIERSTRASS_G2_HPP_
#define SHORT_WEIERSTRASS_G2_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/fields/fp.hpp>

namespace libff {

template<typename SWParamsT>
class short_weierstrass_G2 {
private:
    typename SWParamsT::twist_field X_, Y_, Z_;
public:
    typedef typename SWParamsT::Fq base_field;
    typedef typename SWParamsT::twist_field twist_field;
    typedef typename SWParamsT::Fr scalar_field;

#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static short_weierstrass_G2<SWParamsT> G2_zero;
    static short_weierstrass_G2<SWParamsT> G2_one;
    static twist_field twist;
    static twist_field coeff_a;
    static twist_field coeff_b;
    static void init();

    // using projective coordinates
    short_weierstrass_G2<SWParamsT>() : X_(G2_zero.X_), Y_(G2_zero.Y_), Z_(G2_zero.Z_) {};
    short_weierstrass_G2<SWParamsT>(const twist_field& X, const twist_field& Y, const twist_field& Z) : X_(X), Y_(Y), Z_(Z) {};

    twist_field X() const { return X_; }
    twist_field Y() const { return Y_; }
    twist_field Z() const { return Z_; }

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const short_weierstrass_G2<SWParamsT> &other) const;
    bool operator!=(const short_weierstrass_G2<SWParamsT> &other) const;

    short_weierstrass_G2<SWParamsT> operator+(const short_weierstrass_G2<SWParamsT> &other) const;
    short_weierstrass_G2<SWParamsT> operator-() const;
    short_weierstrass_G2<SWParamsT> operator-(const short_weierstrass_G2<SWParamsT> &other) const;

    short_weierstrass_G2<SWParamsT> add(const short_weierstrass_G2<SWParamsT> &other) const;
    short_weierstrass_G2<SWParamsT> mixed_add(const short_weierstrass_G2<SWParamsT> &other) const;
    short_weierstrass_G2<SWParamsT> dbl() const;
    short_weierstrass_G2<SWParamsT> mul_by_q() const;

    bool is_well_formed() const;

    static short_weierstrass_G2<SWParamsT> zero();
    static short_weierstrass_G2<SWParamsT> one();
    static short_weierstrass_G2<SWParamsT> random_element();

    static size_t size_in_bits() { return twist_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    template<typename SWParamsTT>
    friend std::ostream& operator<<(std::ostream &out, const short_weierstrass_G2<SWParamsTT> &g);
    template<typename SWParamsTT>
    friend std::istream& operator>>(std::istream &in, short_weierstrass_G2<SWParamsTT> &g);

    static void batch_to_special_all_non_zeros(std::vector<short_weierstrass_G2<SWParamsT>> &vec);
};

template<mp_size_t m, typename SWParamsT>
short_weierstrass_G2<SWParamsT> operator*(const bigint<m> &lhs, const short_weierstrass_G2<SWParamsT> &rhs)
{
    return scalar_mul<short_weierstrass_G2<SWParamsT>, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p, typename SWParamsT>
short_weierstrass_G2<SWParamsT> operator*(const Fp_model<m,modulus_p> &lhs, const short_weierstrass_G2<SWParamsT> &rhs)
{
    return scalar_mul<short_weierstrass_G2<SWParamsT>, m>(rhs, lhs.as_bigint());
}

// Begin implementation
template<typename SWParamsT>
using SWtwist_field = typename SWParamsT::twist_field;
template<typename SWParamsT>
using SWFq = typename SWParamsT::Fq;
template<typename SWParamsT>
using SWFr = typename SWParamsT::Fr;

#ifdef PROFILE_OP_COUNTS
template<typename SWParamsT>
long long short_weierstrass_G2<SWParamsT>::add_cnt = 0;
template<typename SWParamsT>
long long short_weierstrass_G2<SWParamsT>::dbl_cnt = 0;
#endif

template<typename SWParamsT>
std::vector<size_t> short_weierstrass_G2<SWParamsT>::wnaf_window_table;
template<typename SWParamsT>
std::vector<size_t> short_weierstrass_G2<SWParamsT>::fixed_base_exp_window_table;
template<typename SWParamsT>
SWtwist_field<SWParamsT> short_weierstrass_G2<SWParamsT>::twist;
template<typename SWParamsT>
SWtwist_field<SWParamsT> short_weierstrass_G2<SWParamsT>::coeff_a;
template<typename SWParamsT>
SWtwist_field<SWParamsT> short_weierstrass_G2<SWParamsT>::coeff_b;
template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::G2_zero;
template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::G2_one;

template<typename SWParamsT>
void short_weierstrass_G2<SWParamsT>::init()
{
    coeff_a = SWParamsT::twist_coeff_a;
    coeff_b = SWParamsT::twist_coeff_b;
    G2_zero = short_weierstrass_G2<SWParamsT>(SWParamsT::G2_zero_X, SWParamsT::G2_zero_Y, SWParamsT::G2_zero_Z);
    G2_one = short_weierstrass_G2<SWParamsT>(SWParamsT::G2_one_X, SWParamsT::G2_one_Y, SWParamsT::G2_one_Z);
}

template<typename SWParamsT>
void short_weierstrass_G2<SWParamsT>::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        short_weierstrass_G2<SWParamsT> copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd*z + %Nd , %Nd*z + %Nd)\n",
                   copy.X_.c1.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   copy.X_.c0.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   copy.Y_.c1.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   copy.Y_.c0.as_bigint().data, SWFq<SWParamsT>::num_limbs);
    }
}

template<typename SWParamsT>
void short_weierstrass_G2<SWParamsT>::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd*z + %Nd : %Nd*z + %Nd : %Nd*z + %Nd)\n",
                   this->X_.c1.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   this->X_.c0.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   this->Y_.c1.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   this->Y_.c0.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   this->Z_.c1.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   this->Z_.c0.as_bigint().data, SWFq<SWParamsT>::num_limbs);
    }
}

template<typename SWParamsT>
void short_weierstrass_G2<SWParamsT>::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X_ = SWtwist_field<SWParamsT>::zero();
        this->Y_ = SWtwist_field<SWParamsT>::one();
        this->Z_ = SWtwist_field<SWParamsT>::zero();
    }
    else
    {
        const SWtwist_field<SWParamsT> Z_inv = Z_.inverse();
        X_ = X_ * Z_inv;
        Y_ = Y_ * Z_inv;
        Z_ = SWtwist_field<SWParamsT>::one();
    }
}

template<typename SWParamsT>
void short_weierstrass_G2<SWParamsT>::to_special()
{
    this->to_affine_coordinates();
}

template<typename SWParamsT>
bool short_weierstrass_G2<SWParamsT>::is_special() const
{
    return (this->is_zero() || this->Z_ == SWtwist_field<SWParamsT>::one());
}

template<typename SWParamsT>
bool short_weierstrass_G2<SWParamsT>::is_zero() const
{
    return (this->X_.is_zero() && this->Z_.is_zero());
}

template<typename SWParamsT>
bool short_weierstrass_G2<SWParamsT>::operator==(const short_weierstrass_G2<SWParamsT> &other) const
{
    if (this->is_zero())
    {
        return other.is_zero();
    }

    if (other.is_zero())
    {
        return false;
    }

    /* now neither is O */

    // X1/Z1 = X2/Z2 <=> X1*Z2 = X2*Z1
    if ((this->X_ * other.Z_) != (other.X_ * this->Z_))
    {
        return false;
    }

    // Y1/Z1 = Y2/Z2 <=> Y1*Z2 = Y2*Z1
    if ((this->Y_ * other.Z_) != (other.Y_ * this->Z_))
    {
        return false;
    }

    return true;
}

template<typename SWParamsT>
bool short_weierstrass_G2<SWParamsT>::operator!=(const short_weierstrass_G2<SWParamsT>& other) const
{
    return !(operator==(other));
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::operator+(const short_weierstrass_G2<SWParamsT> &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return *this;
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // handle double case, and then all the rest
    /*
      The code below is equivalent to (but faster than) the snippet below:

      if (this->operator==(other))
      {
      return this->dbl();
      }
      else
      {
      return this->add(other);
      }
    */

    const SWtwist_field<SWParamsT> X1Z2 = (this->X_) * (other.Z_);        // X1Z2 = X1*Z2
    const SWtwist_field<SWParamsT> X2Z1 = (this->Z_) * (other.X_);        // X2Z1 = X2*Z1

    // (used both in add and double checks)

    const SWtwist_field<SWParamsT> Y1Z2 = (this->Y_) * (other.Z_);        // Y1Z2 = Y1*Z2
    const SWtwist_field<SWParamsT> Y2Z1 = (this->Z_) * (other.Y_);        // Y2Z1 = Y2*Z1

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1)
    {
        // perform dbl case
        const SWtwist_field<SWParamsT> XX   = (this->X_).squared();                   // XX  = X1^2
        const SWtwist_field<SWParamsT> ZZ   = (this->Z_).squared();                   // ZZ  = Z1^2
        const SWtwist_field<SWParamsT> w    = SWParamsT::mul_by_a(ZZ) + (XX + XX + XX); // w   = a*ZZ + 3*XX
        const SWtwist_field<SWParamsT> Y1Z1 = (this->Y_) * (this->Z_);
        const SWtwist_field<SWParamsT> s    = Y1Z1 + Y1Z1;                            // s   = 2*Y1*Z1
        const SWtwist_field<SWParamsT> ss   = s.squared();                            // ss  = s^2
        const SWtwist_field<SWParamsT> sss  = s * ss;                                 // sss = s*ss
        const SWtwist_field<SWParamsT> R    = (this->Y_) * s;                         // R   = Y1*s
        const SWtwist_field<SWParamsT> RR   = R.squared();                            // RR  = R^2
        const SWtwist_field<SWParamsT> B    = ((this->X_)+R).squared()-XX-RR;         // B   = (X1+R)^2 - XX - RR
        const SWtwist_field<SWParamsT> h    = w.squared() - (B+B);                    // h   = w^2 - 2*B
        const SWtwist_field<SWParamsT> X3   = h * s;                                  // X3  = h*s
        const SWtwist_field<SWParamsT> Y3   = w * (B-h)-(RR+RR);                      // Y3  = w*(B-h) - 2*RR
        const SWtwist_field<SWParamsT> Z3   = sss;                                    // Z3  = sss

        return short_weierstrass_G2<SWParamsT>(X3, Y3, Z3);
    }

    // if we have arrived here we are in the add case
    const SWtwist_field<SWParamsT> Z1Z2 = (this->Z_) * (other.Z_);      // Z1Z2 = Z1*Z2
    const SWtwist_field<SWParamsT> u    = Y2Z1 - Y1Z2;                  // u    = Y2*Z1-Y1Z2
    const SWtwist_field<SWParamsT> uu   = u.squared();                  // uu   = u^2
    const SWtwist_field<SWParamsT> v    = X2Z1 - X1Z2;                  // v    = X2*Z1-X1Z2
    const SWtwist_field<SWParamsT> vv   = v.squared();                  // vv   = v^2
    const SWtwist_field<SWParamsT> vvv  = v * vv;                       // vvv  = v*vv
    const SWtwist_field<SWParamsT> R    = vv * X1Z2;                    // R    = vv*X1Z2
    const SWtwist_field<SWParamsT> A    = uu * Z1Z2 - (vvv + R + R);    // A    = uu*Z1Z2 - vvv - 2*R
    const SWtwist_field<SWParamsT> X3   = v * A;                        // X3   = v*A
    const SWtwist_field<SWParamsT> Y3   = u * (R-A) - vvv * Y1Z2;       // Y3   = u*(R-A) - vvv*Y1Z2
    const SWtwist_field<SWParamsT> Z3   = vvv * Z1Z2;                   // Z3   = vvv*Z1Z2

    return short_weierstrass_G2<SWParamsT>(X3, Y3, Z3);
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::operator-() const
{
    return short_weierstrass_G2<SWParamsT>(this->X_, -(this->Y_), this->Z_);
}


template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::operator-(const short_weierstrass_G2<SWParamsT> &other) const
{
    return (*this) + (-other);
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::add(const short_weierstrass_G2<SWParamsT> &other) const
{
    // handle special cases having to do with O
    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return (*this);
    }

    // no need to handle points of order 2,4
    // (they cannot exist in a prime-order subgroup)

    // handle double case
    if (this->operator==(other))
    {
        return this->dbl();
    }

#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2

    const SWtwist_field<SWParamsT> Y1Z2 = (this->Y_) * (other.Z_);        // Y1Z2 = Y1*Z2
    const SWtwist_field<SWParamsT> X1Z2 = (this->X_) * (other.Z_);        // X1Z2 = X1*Z2
    const SWtwist_field<SWParamsT> Z1Z2 = (this->Z_) * (other.Z_);        // Z1Z2 = Z1*Z2
    const SWtwist_field<SWParamsT> u    = (other.Y_) * (this->Z_) - Y1Z2; // u    = Y2*Z1-Y1Z2
    const SWtwist_field<SWParamsT> uu   = u.squared();                    // uu   = u^2
    const SWtwist_field<SWParamsT> v    = (other.X_) * (this->Z_) - X1Z2; // v    = X2*Z1-X1Z2
    const SWtwist_field<SWParamsT> vv   = v.squared();                    // vv   = v^2
    const SWtwist_field<SWParamsT> vvv  = v * vv;                         // vvv  = v*vv
    const SWtwist_field<SWParamsT> R    = vv * X1Z2;                      // R    = vv*X1Z2
    const SWtwist_field<SWParamsT> A    = uu * Z1Z2 - (vvv + R + R);      // A    = uu*Z1Z2 - vvv - 2*R
    const SWtwist_field<SWParamsT> X3   = v * A;                          // X3   = v*A
    const SWtwist_field<SWParamsT> Y3   = u * (R-A) - vvv * Y1Z2;         // Y3   = u*(R-A) - vvv*Y1Z2
    const SWtwist_field<SWParamsT> Z3   = vvv * Z1Z2;                     // Z3   = vvv*Z1Z2

    return short_weierstrass_G2<SWParamsT>(X3, Y3, Z3);
}

template <typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::mixed_add(const short_weierstrass_G2<SWParamsT> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2
    //assert(other.Z == SWtwist_field<SWParamsT>::one());

    if (this->is_zero())
    {
        return other;
    }

    if (other.is_zero())
    {
        return (*this);
    }

#ifdef DEBUG
    assert(other.is_special());
#endif

    const SWtwist_field<SWParamsT> &X1Z2 = (this->X_);                   // X1Z2 = X1*Z2 (but other is special and not zero)
    const SWtwist_field<SWParamsT> X2Z1 = (this->Z_) * (other.X_);       // X2Z1 = X2*Z1

    // (used both in add and double checks)

    const SWtwist_field<SWParamsT> &Y1Z2 = (this->Y_);                   // Y1Z2 = Y1*Z2 (but other is special and not zero)
    const SWtwist_field<SWParamsT> Y2Z1 = (this->Z_) * (other.Y_);       // Y2Z1 = Y2*Z1

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1)
    {
        return this->dbl();
    }

    const SWtwist_field<SWParamsT> u = Y2Z1 - this->Y_;              // u = Y2*Z1-Y1
    const SWtwist_field<SWParamsT> uu = u.squared();                 // uu = u2
    const SWtwist_field<SWParamsT> v = X2Z1 - this->X_;              // v = X2*Z1-X1
    const SWtwist_field<SWParamsT> vv = v.squared();                 // vv = v2
    const SWtwist_field<SWParamsT> vvv = v*vv;                       // vvv = v*vv
    const SWtwist_field<SWParamsT> R = vv * this->X_;                // R = vv*X1
    const SWtwist_field<SWParamsT> A = uu * this->Z_ - vvv - R - R;  // A = uu*Z1-vvv-2*R
    const SWtwist_field<SWParamsT> X3 = v * A;                       // X3 = v*A
    const SWtwist_field<SWParamsT> Y3 = u*(R-A) - vvv * this->Y_;    // Y3 = u*(R-A)-vvv*Y1
    const SWtwist_field<SWParamsT> Z3 = vvv * this->Z_;              // Z3 = vvv*Z1

    return short_weierstrass_G2<SWParamsT>(X3, Y3, Z3);
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::dbl() const
{
#ifdef PROFILE_OP_COUNTS
    this->dbl_cnt++;
#endif
    if (this->is_zero())
    {
        return (*this);
    }
    else
    {
        // NOTE: does not handle O and pts of order 2,4
        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl

        const SWtwist_field<SWParamsT> XX   = (this->X_).squared();                     // XX  = X1^2
        const SWtwist_field<SWParamsT> ZZ   = (this->Z_).squared();                     // ZZ  = Z1^2
        const SWtwist_field<SWParamsT> w    = SWParamsT::mul_by_a(ZZ) + (XX + XX + XX);   // w   = a*ZZ + 3*XX
        const SWtwist_field<SWParamsT> Y1Z1 = (this->Y_) * (this->Z_);
        const SWtwist_field<SWParamsT> s    = Y1Z1 + Y1Z1;                              // s   = 2*Y1*Z1
        const SWtwist_field<SWParamsT> ss   = s.squared();                              // ss  = s^2
        const SWtwist_field<SWParamsT> sss  = s * ss;                                   // sss = s*ss
        const SWtwist_field<SWParamsT> R    = (this->Y_) * s;                           // R   = Y1*s
        const SWtwist_field<SWParamsT> RR   = R.squared();                              // RR  = R^2
        const SWtwist_field<SWParamsT> B    = ((this->X_)+R).squared()-XX-RR;           // B   = (X1+R)^2 - XX - RR
        const SWtwist_field<SWParamsT> h    = w.squared() - (B+B);                      // h   = w^2-2*B
        const SWtwist_field<SWParamsT> X3   = h * s;                                    // X3  = h*s
        const SWtwist_field<SWParamsT> Y3   = w * (B-h)-(RR+RR);                        // Y3  = w*(B-h) - 2*RR
        const SWtwist_field<SWParamsT> Z3   = sss;                                      // Z3  = sss

        return short_weierstrass_G2<SWParamsT>(X3, Y3, Z3);
    }
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::mul_by_q() const
{
    return short_weierstrass_G2<SWParamsT>(
        SWParamsT::twist_mul_by_q_X * (this->X_).Frobenius_map(1),
        SWParamsT::twist_mul_by_q_Y * (this->Y_).Frobenius_map(1),
        (this->Z_).Frobenius_map(1));
}

template<typename SWParamsT>
bool short_weierstrass_G2<SWParamsT>::is_well_formed() const
{
    if (this->is_zero())
    {
        return true;
    }
    else
    {
        /*
          y^2 = x^3 + ax + b

          We are using projective, so equation we need to check is actually

          (y/z)^2 = (x/z)^3 + a (x/z) + b
          z y^2 = x^3  + a z^2 x + b z^3

          z (y^2 - b z^2) = x ( x^2 + a z^2)
        */
        const SWtwist_field<SWParamsT> X2 = this->X_.squared();
        const SWtwist_field<SWParamsT> Y2 = this->Y_.squared();
        const SWtwist_field<SWParamsT> Z2 = this->Z_.squared();
        const SWtwist_field<SWParamsT> aZ2 =  SWParamsT::twist_coeff_a * Z2;

        return (this->Z_ * (Y2 - SWParamsT::twist_coeff_b * Z2) == this->X_ * (X2 + aZ2));
    }
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::zero()
{
    return G2_zero;
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::one()
{
    return G2_one;
}

template<typename SWParamsT>
short_weierstrass_G2<SWParamsT> short_weierstrass_G2<SWParamsT>::random_element()
{
    return (SWFr<SWParamsT>::random_element().as_bigint()) * G2_one;
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream &out, const short_weierstrass_G2<SWParamsT> &g)
{
    short_weierstrass_G2<SWParamsT> copy(g);
    copy.to_affine_coordinates();

    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
#ifdef NO_PT_COMPRESSION
    out << copy.X_ << OUTPUT_SEPARATOR << copy.Y_;
#else
    /* storing LSB of Y */
    out << copy.X_ << OUTPUT_SEPARATOR << (copy.Y_.c0.as_bigint().data[0] & 1);
#endif

    return out;
}

template<typename SWParamsT>
std::istream& operator>>(std::istream &in, short_weierstrass_G2<SWParamsT> &g)
{
    char is_zero;
    SWtwist_field<SWParamsT> tX, tY;

#ifdef NO_PT_COMPRESSION
    in >> is_zero >> tX >> tY;
    is_zero -= '0';
#else
    in.read((char*)&is_zero, 1); // this reads is_zero;
    is_zero -= '0';
    consume_OUTPUT_SEPARATOR(in);

    unsigned char Y_lsb;
    in >> tX;
    consume_OUTPUT_SEPARATOR(in);
    in.read((char*)&Y_lsb, 1);
    Y_lsb -= '0';

    // y = +/- sqrt(x^3 + a*x + b)
    if (!is_zero)
    {
        SWtwist_field<SWParamsT> tX2 = tX.squared();
        SWtwist_field<SWParamsT> tY2 = (tX2 + SWParamsT::twist_coeff_a ) * tX + SWParamsT::twist_coeff_b;
        tY = tY2.sqrt();

        if ((tY.c0.as_bigint().data[0] & 1) != Y_lsb)
        {
            tY = -tY;
        }
    }
#endif
    // using projective coordinates
    if (!is_zero)
    {
        g.X_ = tX;
        g.Y_ = tY;
        g.Z_ = SWtwist_field<SWParamsT>::one();
    }
    else
    {
        g = short_weierstrass_G2<SWParamsT>::zero();
    }

    return in;
}

template<typename SWParamsT>
void short_weierstrass_G2<SWParamsT>::batch_to_special_all_non_zeros(std::vector<short_weierstrass_G2<SWParamsT>> &vec)
{
    std::vector<SWtwist_field<SWParamsT>> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z());
    }
    batch_invert<SWtwist_field<SWParamsT>>(Z_vec);

    const SWtwist_field<SWParamsT> one = SWtwist_field<SWParamsT>::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        vec[i] = short_weierstrass_G2<SWParamsT>(vec[i].X() * Z_vec[i], vec[i].Y() * Z_vec[i], one);
    }
}

} // libff

#endif // SHORT_WEIERSTRASS_G2_HPP_
