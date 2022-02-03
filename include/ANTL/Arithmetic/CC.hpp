/**
 * @file CC.hpp
 * @author Renni deGraaf
 * @version $Header$
 * @brief Template complex number arithmetic class
 */

#ifndef ANTL_CC_H
#define ANTL_CC_H

#include <ANTL/common.hpp>

namespace ANTL
{
  // class prototype needed by friend function definitions
  template < class T > class CC;

  // friend function prototypes needed by friend declatations
  template <class T> void swap(CC<T>& x, CC<T>& y);
  template <class T> bool operator == (const CC<T>& A, const CC<T>& B);
  template <class T> bool operator == (const T& N, const CC<T>& B);
  template <class T> bool operator == (const CC<T>& A, const T& N);
  template <class T> bool operator != (const CC<T>& A, const CC<T>& B);
  template <class T> bool operator != (const T& N, const CC<T>& B);
  template <class T> bool operator != (const CC<T>& A, const T& N);
  template <class T> void add(CC<T>& C, const CC<T>& A, const CC<T>& B);
  template <class T> void add(CC<T>& C, const CC<T>& A, const T& B);
  template <class T> void add(CC<T>& C, const T& A, const T& B);
  template <class T> void subtract(CC<T>& C, const CC<T>& A, const CC<T>& B);
  template <class T> void subtract(CC<T>& C, const CC<T>& A, const T& B);
  template <class T> void subtract(CC<T>& C, const T& A, const T& B);
  template <class T> void multiply(CC<T>& C, const CC<T>& A, const CC<T>& B);
  template <class T> void multiply(CC<T>& C, const CC<T>& A, const T& B);
  template <class T> void multiply(CC<T>& C, const T& A, const CC<T>& B);
  template <class T> void multiply(CC<T>& C, const T& A, const T& B);
  template <class T> void divide(CC<T>& C, const CC<T>& A, const CC<T>& B);
  template <class T> void divide(CC<T>& C, const CC<T>& A, const T& B);
  template <class T> void divide(CC<T>& C, const T& A, const CC<T>& B);
  template <class T> void divide(CC<T>& C, const T& A, const T& B);
  template <class T> void square(CC<T>& C, const CC<T>& A);
  template <class T> void square(CC<T>& C, const T& A);
  template <class T> void power(CC<T>& C, const CC<T>& A, const unsigned long e);
  template <class T> void power(CC<T>& C, const CC<T>& A, const ZZ& e);
  template <class T> CC<T> operator - (const CC<T>& A);
  template <class T> CC<T> operator + (const CC<T>& A, const CC<T>& B);
  template <class T> CC<T> operator + (const CC<T>& A, const T& B);
  template <class T> CC<T> operator - (const CC<T>& A, const CC<T>& B);
  template <class T> CC<T> operator - (const CC<T>& A, const T& B);
  template <class T> CC<T> operator * (const CC<T>& A, const CC<T>& B);
  template <class T> CC<T> operator * (const CC<T>& A, const T& B);
  template <class T> CC<T> operator * (const T& A, const CC<T>& B);
  template <class T> CC<T> operator / (const CC<T>& A, const CC<T>& B);
  template <class T> CC<T> operator / (const CC<T>& A, const T& B);
  template <class T> CC<T> operator / (const T& A, const CC<T>& B);
  template <class T> const CC<T>& operator += (CC<T>& A, const CC<T>& B);
  template <class T> const CC<T>& operator -= (CC<T>& A, const CC<T>& B);
  template <class T> const CC<T>& operator *= (CC<T>& A, const CC<T>& B);
  template <class T> const CC<T>& operator /= (CC<T>& A, const CC<T>& B);
  template <class T> std::istream& operator >> (std::istream& in, CC<T>& X);
  template <class T> void inverse(CC<T>& C, const CC<T>& A);
  template <class T> void conjugate(CC<T>& C, const CC<T>& A);
  template <class T> std::ostream& operator << (std::ostream& out, const CC<T>& X);

  template <class T>
  class CC
  {
  protected:
    T r;				// real part
    T i;				// imaginary part

  public:
    CC ();
    CC(const CC<T>& x);
    CC(const T& x);
    CC(const T& real, const T& imag);
    ~CC ();

    /* access */
    void get(T& real, T& imag);
    const T& real() const;                        // not in quadratic_number
    const T& get_r() const;
    const T& imaginary() const;                   // not in quadratic_number
    const T& get_i() const;

    /* assignment */
    void assign (const T& newreal, const T& newimag);
    void assign(const T& newreal);
    void assign(const CC<T>& x);
    void clear (void);                            // not in quadratic_number
    void assign_zero();
    void set (void);                              // not in quadratic_number
    void assign_one();
    void set_i (void);                            // not in quadratic_number
    friend void swap<T>(CC<T>& x, CC<T>& y);
    const CC<T>& operator = (const CC<T>& A);
    const CC<T>& operator = (const T& A);         // not in quadratic_number

    /* comparisions */
    bool IsZero();                                // not in quadratic_number
    bool IsOne();                                 // not in quadratic_number
    bool is_zero();
    bool is_one();
    bool is_i();                                  // not in quadratic_number
    bool is_real (void);                          // not in quadratic_number
    bool is_equal(const CC<T>& B);
    bool is_equal(const T& ea);
    bool is_equal(const T& ea, const T& eb);
    friend bool operator == <T> (const CC<T>& A, const CC<T>& B);
    friend bool operator == <T> (const T& N, const CC<T>& B);
    friend bool operator == <T> (const CC<T>& A, const T& N);
    friend bool operator != <T> (const CC<T>& A, const CC<T>& B);
    friend bool operator != <T> (const T& N, const CC<T>& B);
    friend bool operator != <T> (const CC<T>& A, const T& N);

    /* arithmetic */
    T norm() const;                               // not in quadratic_number
    T trace() const;
    void negate();
    void negate(CC<T>& B);                        // not in quadratic_number
    void invert();
    void invert(CC<T>& B);                        // not in quadratic_number
    friend void inverse <T> (CC<T>& C, const CC<T>& A);
    friend void conjugate <T> (CC<T>& C, const CC<T>& A);

    friend void add <T> (CC<T>& C, const CC<T>& A, const CC<T>& B);
    friend void add <T> (CC<T>& C, const CC<T>& A, const T& B); // not in quadratic_number
    friend void add <T> (CC<T>& C, const T& A, const T& B); // not in quadratic_number

    friend void subtract <T> (CC<T>& C, const CC<T>& A, const CC<T>& B);
    friend void subtract <T> (CC<T>& C, const CC<T>& A, const T& B); // not in quadratic_number
    friend void subtract <T> (CC<T>& C, const T& A, const T& B); // not in quadratic_number

    friend void multiply <T> (CC<T>& C, const CC<T>& A, const CC<T>& B);
    friend void multiply <T> (CC<T>& C, const CC<T>& A, const T& B);
    friend void multiply <T> (CC<T>& C, const T& A, const CC<T>& B);
    friend void multiply <T> (CC<T>& C, const T& A, const T& B); // not in quadratic_number

    friend void divide <T> (CC<T>& C, const CC<T>& A, const CC<T>& B);
    friend void divide <T> (CC<T>& C, const CC<T>& A, const T& B);
    friend void divide <T> (CC<T>& C, const T& A, const CC<T>& B);
    friend void divide <T> (CC<T>& C, const T& A, const T& B); // not in quadratic_number

    friend void square <T> (CC<T>& C, const CC<T>& A);
    friend void square <T> (CC<T>& C, const T& A); // not in quadratic_number

    friend void power <T> (CC<T>& C, const CC<T>& A, const unsigned long e);
    friend void power <T> (CC<T>& C, const CC<T>& A, const ZZ& e);

    friend CC<T> operator - <T> (const CC<T>& A);
    friend CC<T> operator + <T> (const CC<T>& A, const CC<T>& B);
    friend CC<T> operator + <T> (const CC<T>& A, const T& B); // not in quadratic_number
    friend CC<T> operator - <T> (const CC<T>& A, const CC<T>& B);
    friend CC<T> operator - <T> (const CC<T>& A, const T& B); // not in quadratic_number
    friend CC<T> operator * <T> (const CC<T>& A, const CC<T>& B);
    friend CC<T> operator * <T> (const CC<T>& A, const T& B);
    friend CC<T> operator * <T> (const T& A, const CC<T>& B);
    friend CC<T> operator / <T> (const CC<T>& A, const CC<T>& B);
    friend CC<T> operator / <T> (const CC<T>& A, const T& B);
    friend CC<T> operator / <T> (const T& A, const CC<T>& B);
    friend const CC<T>& operator += <T> (CC<T>& A, const CC<T>& B);
    friend const CC<T>& operator -= <T> (CC<T>& A, const CC<T>& B);
    friend const CC<T>& operator *= <T> (CC<T>& A, const CC<T>& B);
    friend const CC<T>& operator /= <T> (CC<T>& A, const CC<T>& B);

    /* input/output */
    friend std::istream& operator >> <T> (std::istream& in, CC<T>& X);
    friend std::ostream& operator << <T> (std::ostream& out, const CC<T>& X);

    /* precision */
    void SetPrecision (long p);                   // not in quadratic_number
    long Precision (void);                        // not in quadratic_number
    void SetOutputPrecision (long p);             // not in quadratic_number
    long OutputPrecision (void);                  // not in quadratic_number
  };

  //
  // Template-friendly conversion methods.
  // Note that xdouble and quad_float are incompatible.
  //
  template <class T> T to(const CC<float>& a) { return T(a); }
  template <class T> T to(const CC<double>& a) { return T(a); }
  //template <class T> T to(const CC<xdouble>& a) { return T(a); }
  template <class T> T to(const CC<quad_float>& a) { return T(a); }
  template <class T> T to(const CC<RR>& a) { return T(a); }

  template<> inline CC<float> to< CC<float> > (const CC<float> & a)
  { return a; }
  template<> inline CC<float> to< CC<float> > (const CC<double> & a)
  { return CC<float>(to<float>(a.real()), to<float>(a.imaginary())); }
  //template<> inline CC<float> to< CC<float> > (const CC<xdouble> & a)
  //{ return CC<float>(to<float>(a.real()), to<float>(a.imaginary())); }
  template<> inline CC<float> to< CC<float> > (const CC<quad_float> & a)
  { return CC<float>(to<float>(a.real()), to<float>(a.imaginary())); }
  template<> inline CC<float> to< CC<float> > (const CC<RR> & a)
  { return CC<float>(to<float>(a.real()), to<float>(a.imaginary())); }

  template<> inline CC<double> to< CC<double> > (const CC<float> & a)
  { return CC<double>(to<double>(a.real()), to<double>(a.imaginary())); }
  template<> inline CC<double> to< CC<double> > (const CC<double> & a)
  { return a; }
  //template<> inline CC<double> to< CC<double> > (const CC<xdouble> & a)
  //{ return CC<double>(to<double>(a.real()), to<double>(a.imaginary())); }
  template<> inline CC<double> to< CC<double> > (const CC<quad_float> & a)
  { return CC<double>(to<double>(a.real()), to<double>(a.imaginary())); }
  template<> inline CC<double> to< CC<double> > (const CC<RR> & a)
  { return CC<double>(to<double>(a.real()), to<double>(a.imaginary())); }

  template<> inline CC<xdouble> to< CC<xdouble> > (const CC<float> & a)
  { return CC<xdouble>(to<xdouble>(a.real()), to<xdouble>(a.imaginary())); }
  template<> inline CC<xdouble> to< CC<xdouble> > (const CC<double> & a)
  { return CC<xdouble>(to<xdouble>(a.real()), to<xdouble>(a.imaginary())); }
  //template<> inline CC<xdouble> to< CC<xdouble> > (const CC<xdouble> & a)
  //{ return a; }
  //template<> inline CC<xdouble> to< CC<xdouble> > (const CC<RR> & a)
  //{ return CC<xdouble>(to<xdouble>(a.real()), to<xdouble>(a.imaginary())); }

  template<> inline CC<quad_float> to< CC<quad_float> > (const CC<float> & a)
  { return CC<quad_float>(to<quad_float>(a.real()), to<quad_float>(a.imaginary())); }
  template<> inline CC<quad_float> to< CC<quad_float> > (const CC<double> & a)
  { return CC<quad_float>(to<quad_float>(a.real()), to<quad_float>(a.imaginary())); }
  template<> inline CC<quad_float> to< CC<quad_float> > (const CC<quad_float> & a)
  { return a; }
  template<> inline CC<quad_float> to< CC<quad_float> > (const CC<RR> & a)
  { return CC<quad_float>(to<quad_float>(a.real()), to<quad_float>(a.imaginary())); }

  template<> inline CC<RR> to< CC<RR> > (const CC<float> & a)
  { return CC<RR>(to<RR>(a.real()), to<RR>(a.imaginary())); }
  template<> inline CC<RR> to< CC<RR> > (const CC<double> & a)
  { return CC<RR>(to<RR>(a.real()), to<RR>(a.imaginary())); }
  //template<> inline CC<RR> to< CC<RR> > (const CC<xdouble> & a)
  //{ return CC<RR>(to<RR>(a.real()), to<RR>(a.imaginary())); }
  template<> inline CC<RR> to< CC<RR> > (const CC<quad_float> & a)
  { return CC<RR>(to<RR>(a.real()), to<RR>(a.imaginary())); }
  template<> inline CC<RR> to< CC<RR> > (const CC<RR> & a)
  { return a; }

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Arithmetic/CC_impl.hpp"

#endif // guard
