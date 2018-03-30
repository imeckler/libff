/** @file
 *****************************************************************************

 Declaration of interfaces for the MNT4 G1 group.

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef SHORT_WEIERSTRASS_G1_HPP_
#define SHORT_WEIERSTRASS_G1_HPP_

#include <vector>

#include <libff/algebra/curves/curve_utils.hpp>
#include <libff/algebra/fields/fp.hpp>

namespace libff {

template<typename SWParamsT>
using SWFq = typename SWParamsT::Fq;

template<typename SWParamsT>
using SWFr = typename SWParamsT::Fr;

template<typename SWParamsT>
class short_weierstrass_G1;

template<typename SWParamsT>
std::ostream& operator<<(std::ostream &, const short_weierstrass_G1<SWParamsT>&);
template<typename SWParamsT>
std::istream& operator>>(std::istream &, short_weierstrass_G1<SWParamsT>&);

template<typename SWParamsT>
class short_weierstrass_G1 {
    typedef typename SWParamsT::Fq base_field;
    typedef typename SWParamsT::Fr scalar_field;

private:
    base_field X_, Y_, Z_;
public:
#ifdef PROFILE_OP_COUNTS
    static long long add_cnt;
    static long long dbl_cnt;
#endif
    static std::vector<size_t> wnaf_window_table;
    static std::vector<size_t> fixed_base_exp_window_table;
    static short_weierstrass_G1<SWParamsT> G1_zero;
    static short_weierstrass_G1<SWParamsT> G1_one;
    static base_field coeff_a;
    static base_field coeff_b;
    static void init();

    // using projective coordinates
    short_weierstrass_G1<SWParamsT>() : X_(SWParamsT::G1_zero_X), Y_(SWParamsT::G1_zero_Y), Z_(SWParamsT::G1_zero_Z) {};
    short_weierstrass_G1<SWParamsT>(const base_field& X, const base_field& Y) : X_(X), Y_(Y), Z_(base_field::one()) {}
    short_weierstrass_G1<SWParamsT>(const base_field& X, const base_field& Y, const base_field& Z) : X_(X), Y_(Y), Z_(Z) {}

    base_field X() const { return X_; }
    base_field Y() const { return Y_; }
    base_field Z() const { return Z_; }

    void print() const;
    void print_coordinates() const;

    void to_affine_coordinates();
    void to_special();
    bool is_special() const;

    bool is_zero() const;

    bool operator==(const short_weierstrass_G1<SWParamsT> &other) const;
    bool operator!=(const short_weierstrass_G1<SWParamsT> &other) const;

    short_weierstrass_G1<SWParamsT> operator+(const short_weierstrass_G1<SWParamsT> &other) const;
    short_weierstrass_G1<SWParamsT> operator-() const;
    short_weierstrass_G1<SWParamsT> operator-(const short_weierstrass_G1<SWParamsT> &other) const;

    short_weierstrass_G1<SWParamsT> add(const short_weierstrass_G1<SWParamsT> &other) const;
    short_weierstrass_G1<SWParamsT> mixed_add(const short_weierstrass_G1<SWParamsT> &other) const;
    short_weierstrass_G1<SWParamsT> dbl() const;

    bool is_well_formed() const;

    static short_weierstrass_G1<SWParamsT> zero();
    static short_weierstrass_G1<SWParamsT> one();
    static short_weierstrass_G1<SWParamsT> random_element();

    static size_t size_in_bits() { return base_field::size_in_bits() + 1; }
    static bigint<base_field::num_limbs> base_field_char() { return base_field::field_char(); }
    static bigint<scalar_field::num_limbs> order() { return scalar_field::field_char(); }

    template <typename SWParamsTT>
    friend std::ostream& operator<<(std::ostream &out, const short_weierstrass_G1<SWParamsTT> &g);
    template <typename SWParamsTT>
    friend std::istream& operator>>(std::istream &in, short_weierstrass_G1<SWParamsTT> &g);

    static void batch_to_special_all_non_zeros(std::vector<short_weierstrass_G1<SWParamsT>> &vec);
};

template<mp_size_t m, typename SWParamsT>
short_weierstrass_G1<SWParamsT> operator*(const bigint<m> &lhs, const short_weierstrass_G1<SWParamsT> &rhs)
{
    return scalar_mul<short_weierstrass_G1<SWParamsT>, m>(rhs, lhs);
}

template<mp_size_t m, const bigint<m>& modulus_p, typename SWParamsT>
short_weierstrass_G1<SWParamsT> operator*(const Fp_model<m,modulus_p> &lhs, const short_weierstrass_G1<SWParamsT> &rhs)
{
    return scalar_mul<short_weierstrass_G1<SWParamsT>, m>(rhs, lhs.as_bigint());
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream& out, const std::vector<short_weierstrass_G1<SWParamsT>> &v);
template<typename SWParamsT>
std::istream& operator>>(std::istream& in, std::vector<short_weierstrass_G1<SWParamsT>> &v);

// Begin implementation

#ifdef PROFILE_OP_COUNTS
template<typename SWParamsT>
long long short_weierstrass_G1<SWParamsT>::add_cnt = 0;
template<typename SWParamsT>
long long short_weierstrass_G1<SWParamsT>::dbl_cnt = 0;
#endif

template<typename SWParamsT>
std::vector<size_t> short_weierstrass_G1<SWParamsT>::wnaf_window_table;
template<typename SWParamsT>
std::vector<size_t> short_weierstrass_G1<SWParamsT>::fixed_base_exp_window_table;

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::G1_zero;
template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::G1_one;
template<typename SWParamsT>
SWFq<SWParamsT> short_weierstrass_G1<SWParamsT>::coeff_a;
template<typename SWParamsT>
SWFq<SWParamsT> short_weierstrass_G1<SWParamsT>::coeff_b;

template<typename SWParamsT>
void short_weierstrass_G1<SWParamsT>::init()
{
    coeff_a = SWParamsT::coeff_a;
    coeff_b = SWParamsT::coeff_b;
    G1_zero = short_weierstrass_G1<SWParamsT>(SWParamsT::G1_zero_X, SWParamsT::G1_zero_Y, SWParamsT::G1_zero_Z);
    G1_one = short_weierstrass_G1<SWParamsT>(SWParamsT::G1_one_X, SWParamsT::G1_one_Y, SWParamsT::G1_one_Z);
}

template<typename SWParamsT>
void short_weierstrass_G1<SWParamsT>::print() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        short_weierstrass_G1<SWParamsT> copy(*this);
        copy.to_affine_coordinates();
        gmp_printf("(%Nd , %Nd)\n",
                   copy.X_.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   copy.Y_.as_bigint().data, SWFq<SWParamsT>::num_limbs);
    }
}

template<typename SWParamsT>
void short_weierstrass_G1<SWParamsT>::print_coordinates() const
{
    if (this->is_zero())
    {
        printf("O\n");
    }
    else
    {
        gmp_printf("(%Nd : %Nd : %Nd)\n",
                   this->X_.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   this->Y_.as_bigint().data, SWFq<SWParamsT>::num_limbs,
                   this->Z_.as_bigint().data, SWFq<SWParamsT>::num_limbs);
    }
}

template<typename SWParamsT>
void short_weierstrass_G1<SWParamsT>::to_affine_coordinates()
{
    if (this->is_zero())
    {
        this->X_ = SWFq<SWParamsT>::zero();
        this->Y_ = SWFq<SWParamsT>::one();
        this->Z_ = SWFq<SWParamsT>::zero();
    }
    else
    {
        const SWFq<SWParamsT> Z_inv = Z_.inverse();
        this->X_ = this->X_ * Z_inv;
        this->Y_ = this->Y_ * Z_inv;
        this->Z_ = SWFq<SWParamsT>::one();
    }
}

template<typename SWParamsT>
void short_weierstrass_G1<SWParamsT>::to_special()
{
    this->to_affine_coordinates();
}

template<typename SWParamsT>
bool short_weierstrass_G1<SWParamsT>::is_special() const
{
    return (this->is_zero() || this->Z_ == SWFq<SWParamsT>::one());
}

template<typename SWParamsT>
bool short_weierstrass_G1<SWParamsT>::is_zero() const
{
    return (this->X_.is_zero() && this->Z_.is_zero());
}

template<typename SWParamsT>
bool short_weierstrass_G1<SWParamsT>::operator==(const short_weierstrass_G1<SWParamsT> &other) const
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
bool short_weierstrass_G1<SWParamsT>::operator!=(const short_weierstrass_G1<SWParamsT>& other) const
{
    return !(operator==(other));
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::operator+(const short_weierstrass_G1<SWParamsT> &other) const
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

    const SWFq<SWParamsT> X1Z2 = (this->X_) * (other.Z_);        // X1Z2 = X1*Z2
    const SWFq<SWParamsT> X2Z1 = (this->Z_) * (other.X_);        // X2Z1 = X2*Z1

    // (used both in add and double checks)

    const SWFq<SWParamsT> Y1Z2 = (this->Y_) * (other.Z_);        // Y1Z2 = Y1*Z2
    const SWFq<SWParamsT> Y2Z1 = (this->Z_) * (other.Y_);        // Y2Z1 = Y2*Z1

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1)
    {
        // perform dbl case
        const SWFq<SWParamsT> XX   = (this->X_).squared();                   // XX  = X1^2
        const SWFq<SWParamsT> ZZ   = (this->Z_).squared();                   // ZZ  = Z1^2
        const SWFq<SWParamsT> w    = SWParamsT::coeff_a * ZZ + (XX + XX + XX); // w   = a*ZZ + 3*XX
        const SWFq<SWParamsT> Y1Z1 = (this->Y_) * (this->Z_);
        const SWFq<SWParamsT> s    = Y1Z1 + Y1Z1;                            // s   = 2*Y1*Z1
        const SWFq<SWParamsT> ss   = s.squared();                            // ss  = s^2
        const SWFq<SWParamsT> sss  = s * ss;                                 // sss = s*ss
        const SWFq<SWParamsT> R    = (this->Y_) * s;                         // R   = Y1*s
        const SWFq<SWParamsT> RR   = R.squared();                            // RR  = R^2
        const SWFq<SWParamsT> B    = ((this->X_)+R).squared()-XX-RR;         // B   = (X1+R)^2 - XX - RR
        const SWFq<SWParamsT> h    = w.squared() - (B+B);                    // h   = w^2 - 2*B
        const SWFq<SWParamsT> X3   = h * s;                                  // X3  = h*s
        const SWFq<SWParamsT> Y3   = w * (B-h)-(RR+RR);                      // Y3  = w*(B-h) - 2*RR
        const SWFq<SWParamsT> Z3   = sss;                                    // Z3  = sss

        return short_weierstrass_G1<SWParamsT>(X3, Y3, Z3);
    }

    // if we have arrived here we are in the add case
    const SWFq<SWParamsT> Z1Z2 = (this->Z_) * (other.Z_);        // Z1Z2 = Z1*Z2
    const SWFq<SWParamsT> u    = Y2Z1 - Y1Z2; // u    = Y2*Z1-Y1Z2
    const SWFq<SWParamsT> uu   = u.squared();                  // uu   = u^2
    const SWFq<SWParamsT> v    = X2Z1 - X1Z2; // v    = X2*Z1-X1Z2
    const SWFq<SWParamsT> vv   = v.squared();                  // vv   = v^2
    const SWFq<SWParamsT> vvv  = v * vv;                       // vvv  = v*vv
    const SWFq<SWParamsT> R    = vv * X1Z2;                    // R    = vv*X1Z2
    const SWFq<SWParamsT> A    = uu * Z1Z2 - (vvv + R + R);    // A    = uu*Z1Z2 - vvv - 2*R
    const SWFq<SWParamsT> X3   = v * A;                        // X3   = v*A
    const SWFq<SWParamsT> Y3   = u * (R-A) - vvv * Y1Z2;       // Y3   = u*(R-A) - vvv*Y1Z2
    const SWFq<SWParamsT> Z3   = vvv * Z1Z2;                   // Z3   = vvv*Z1Z2

    return short_weierstrass_G1<SWParamsT>(X3, Y3, Z3);
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::operator-() const
{
    return short_weierstrass_G1<SWParamsT>(this->X_, -(this->Y_), this->Z_);
}


template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::operator-(const short_weierstrass_G1<SWParamsT> &other) const
{
    return (*this) + (-other);
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::add(const short_weierstrass_G1<SWParamsT> &other) const
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

    const SWFq<SWParamsT> Y1Z2 = (this->Y_) * (other.Z_);        // Y1Z2 = Y1*Z2
    const SWFq<SWParamsT> X1Z2 = (this->X_) * (other.Z_);        // X1Z2 = X1*Z2
    const SWFq<SWParamsT> Z1Z2 = (this->Z_) * (other.Z_);        // Z1Z2 = Z1*Z2
    const SWFq<SWParamsT> u    = (other.Y_) * (this->Z_) - Y1Z2; // u    = Y2*Z1-Y1Z2
    const SWFq<SWParamsT> uu   = u.squared();                    // uu   = u^2
    const SWFq<SWParamsT> v    = (other.X_) * (this->Z_) - X1Z2; // v    = X2*Z1-X1Z2
    const SWFq<SWParamsT> vv   = v.squared();                    // vv   = v^2
    const SWFq<SWParamsT> vvv  = v * vv;                         // vvv  = v*vv
    const SWFq<SWParamsT> R    = vv * X1Z2;                      // R    = vv*X1Z2
    const SWFq<SWParamsT> A    = uu * Z1Z2 - (vvv + R + R);      // A    = uu*Z1Z2 - vvv - 2*R
    const SWFq<SWParamsT> X3   = v * A;                          // X3   = v*A
    const SWFq<SWParamsT> Y3   = u * (R-A) - vvv * Y1Z2;         // Y3   = u*(R-A) - vvv*Y1Z2
    const SWFq<SWParamsT> Z3   = vvv * Z1Z2;                     // Z3   = vvv*Z1Z2

    return short_weierstrass_G1<SWParamsT>(X3, Y3, Z3);
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::mixed_add(const short_weierstrass_G1<SWParamsT> &other) const
{
#ifdef PROFILE_OP_COUNTS
    this->add_cnt++;
#endif
    // NOTE: does not handle O and pts of order 2,4
    // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-add-1998-cmo-2
    //assert(other.Z == mnt4_Fq::one());

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

    const SWFq<SWParamsT> &X1Z2 = (this->X_);                    // X1Z2 = X1*Z2 (but other is special and not zero)
    const SWFq<SWParamsT> X2Z1 = (this->Z_) * (other.X_);        // X2Z1 = X2*Z1

    // (used both in add and double checks)

    const SWFq<SWParamsT> &Y1Z2 = (this->Y_);                    // Y1Z2 = Y1*Z2 (but other is special and not zero)
    const SWFq<SWParamsT> Y2Z1 = (this->Z_) * (other.Y_);        // Y2Z1 = Y2*Z1

    if (X1Z2 == X2Z1 && Y1Z2 == Y2Z1)
    {
        return this->dbl();
    }

    const SWFq<SWParamsT> u = Y2Z1 - this->Y_;              // u = Y2*Z1-Y1
    const SWFq<SWParamsT> uu = u.squared();                 // uu = u2
    const SWFq<SWParamsT> v = X2Z1 - this->X_;              // v = X2*Z1-X1
    const SWFq<SWParamsT> vv = v.squared();                 // vv = v2
    const SWFq<SWParamsT> vvv = v*vv;                       // vvv = v*vv
    const SWFq<SWParamsT> R = vv * this->X_;                // R = vv*X1
    const SWFq<SWParamsT> A = uu * this->Z_ - vvv - R - R;  // A = uu*Z1-vvv-2*R
    const SWFq<SWParamsT> X3 = v * A;                       // X3 = v*A
    const SWFq<SWParamsT> Y3 = u*(R-A) - vvv * this->Y_;    // Y3 = u*(R-A)-vvv*Y1
    const SWFq<SWParamsT> Z3 = vvv * this->Z_;              // Z3 = vvv*Z1

    return short_weierstrass_G1<SWParamsT>(X3, Y3, Z3);
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::dbl() const
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

        const SWFq<SWParamsT> XX   = (this->X_).squared();                   // XX  = X1^2
        const SWFq<SWParamsT> ZZ   = (this->Z_).squared();                   // ZZ  = Z1^2
        const SWFq<SWParamsT> w    = SWParamsT::coeff_a * ZZ + (XX + XX + XX); // w   = a*ZZ + 3*XX
        const SWFq<SWParamsT> Y1Z1 = (this->Y_) * (this->Z_);
        const SWFq<SWParamsT> s    = Y1Z1 + Y1Z1;                            // s   = 2*Y1*Z1
        const SWFq<SWParamsT> ss   = s.squared();                            // ss  = s^2
        const SWFq<SWParamsT> sss  = s * ss;                                 // sss = s*ss
        const SWFq<SWParamsT> R    = (this->Y_) * s;                         // R   = Y1*s
        const SWFq<SWParamsT> RR   = R.squared();                            // RR  = R^2
        const SWFq<SWParamsT> B    = ((this->X_)+R).squared()-XX-RR;         // B   = (X1+R)^2 - XX - RR
        const SWFq<SWParamsT> h    = w.squared() - (B+B);                    // h   = w^2 - 2*B
        const SWFq<SWParamsT> X3   = h * s;                                  // X3  = h*s
        const SWFq<SWParamsT> Y3   = w * (B-h)-(RR+RR);                      // Y3  = w*(B-h) - 2*RR
        const SWFq<SWParamsT> Z3   = sss;                                    // Z3  = sss

        return short_weierstrass_G1<SWParamsT>(X3, Y3, Z3);
    }
}

template<typename SWParamsT>
bool short_weierstrass_G1<SWParamsT>::is_well_formed() const
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
        const SWFq<SWParamsT> X2 = this->X_.squared();
        const SWFq<SWParamsT> Y2 = this->Y_.squared();
        const SWFq<SWParamsT> Z2 = this->Z_.squared();

        return (this->Z_ * (Y2 - SWParamsT::coeff_b * Z2) == this->X_ * (X2 + SWParamsT::coeff_a * Z2));
    }
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::zero()
{
    return short_weierstrass_G1<SWParamsT>::G1_zero;
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::one()
{
    return short_weierstrass_G1<SWParamsT>::G1_one;
}

template<typename SWParamsT>
short_weierstrass_G1<SWParamsT> short_weierstrass_G1<SWParamsT>::random_element()
{
    return (scalar_field::random_element().as_bigint()) * short_weierstrass_G1<SWParamsT>::one();
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream &out, const short_weierstrass_G1<SWParamsT> &g)
{
    short_weierstrass_G1<SWParamsT> copy(g);
    copy.to_affine_coordinates();

    out << (copy.is_zero() ? 1 : 0) << OUTPUT_SEPARATOR;
#ifdef NO_PT_COMPRESSION
    out << copy.X_ << OUTPUT_SEPARATOR << copy.Y_;
#else
    /* storing LSB of Y */
    out << copy.X_ << OUTPUT_SEPARATOR << (copy.Y_.as_bigint().data[0] & 1);
#endif

    return out;
}

template<typename SWParamsT>
std::istream& operator>>(std::istream &in, short_weierstrass_G1<SWParamsT> &g)
{
    char is_zero;
    SWFq<SWParamsT> tX, tY;

#ifdef NO_PT_COMPRESSION
    in >> is_zero >> tX >> tY;
    is_zero -= '0';
#else
    in.read((char*)&is_zero, 1);
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
        SWFq<SWParamsT> tX2 = tX.squared();
        SWFq<SWParamsT> tY2 = (tX2 + SWParamsT::coeff_a) * tX + SWParamsT::coeff_b;
        tY = tY2.sqrt();

        if ((tY.as_bigint().data[0] & 1) != Y_lsb)
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
        g.Z_ = SWFq<SWParamsT>::one();
    }
    else
    {
        g = short_weierstrass_G1<SWParamsT>::zero();
    }

    return in;
}

template<typename SWParamsT>
std::ostream& operator<<(std::ostream& out, const std::vector<short_weierstrass_G1<SWParamsT>> &v)
{
    out << v.size() << "\n";
    for (const short_weierstrass_G1<SWParamsT>& t : v)
    {
        out << t << OUTPUT_NEWLINE;
    }

    return out;
}

template<typename SWParamsT>
std::istream& operator>>(std::istream& in, std::vector<short_weierstrass_G1<SWParamsT>> &v)
{
    v.clear();

    size_t s;
    in >> s;

    consume_newline(in);

    v.reserve(s);

    for (size_t i = 0; i < s; ++i)
    {
        short_weierstrass_G1<SWParamsT> g;
        in >> g;
        consume_OUTPUT_NEWLINE(in);
        v.emplace_back(g);
    }

    return in;
}

template<typename SWParamsT>
void short_weierstrass_G1<SWParamsT>::batch_to_special_all_non_zeros(std::vector<short_weierstrass_G1<SWParamsT>> &vec)
{
    std::vector<SWFq<SWParamsT>> Z_vec;
    Z_vec.reserve(vec.size());

    for (auto &el: vec)
    {
        Z_vec.emplace_back(el.Z());
    }
    batch_invert<SWFq<SWParamsT>>(Z_vec);

    const SWFq<SWParamsT> one = SWFq<SWParamsT>::one();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        vec[i] = short_weierstrass_G1<SWParamsT>(vec[i].X() * Z_vec[i], vec[i].Y() * Z_vec[i], one);
    }
}

} // libff

#endif // SHORT_WEIERSTRASS_G1_HPP_
