/**
 * @file QuadraticIdeal.hpp
 * @version $Header: /scrinium/ANTL/ANTL/include/ANTL/quadratic/QuadraticIdeal.hpp,v 1.1 2013/02/26 06:29:10 jacobs Exp $
 *
 * Basically a quadratic form, but with additional q, so that the ideal is
 * I = (num/dem)(aZ + ((b+Sqrt(Delta))/2) Z)
 */

#ifndef ANTL_QUADRATIC_QUADRATIC_IDEAL_H
#define ANTL_QUADRATIC_QUADRATIC_IDEAL_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Arithmetic/QQ.hpp>

namespace ANTL
{

  template < class T> class QuadraticIdeal;

  template <class T>
  void multiply(QuadraticIdeal<T> & C,
		const QuadraticIdeal<T> & A, const QuadraticIdeal<T> &B);

  template < class T >
  void multiply_real(QuadraticIdeal<T> & C, const QuadraticIdeal<T> & A,
		     const QuadraticIdeal<T> & B, QuadraticNumber<T> & gamma,
		     const QuadraticNumber<T> & alpha, const QuadraticNumber<T> & beta);

  //similarly, needed for squaring, powers
  template <class T>
  void square(QuadraticIdeal<T> & C, const QuadraticIdeal<T> & A);

  template <class T>
  void power (QuadraticIdeal<T> &C, const QuadraticIdeal<T> & A, const ZZ & n);


  template <class T>
  QuadraticIdeal<T> inverse(const QuadraticIdeal<T> & A);

  template <class T>
  QuadraticIdeal<T> operator - (const QuadraticIdeal<T> & A);

  template <class T>
  QuadraticIdeal<T> operator * (const QuadraticIdeal<T> & A,
				 const QuadraticIdeal<T> & B);

  template <class T>
  void swap (QuadraticIdeal<T> & A, QuadraticIdeal<T> & B);

  template <class T>
  bool operator == (const QuadraticIdeal<T> & A, const QuadraticIdeal<T> & B);

  template <class T>
  bool operator != (const QuadraticIdeal<T> & A, const QuadraticIdeal<T> & B);

  template <class T>
  std::istream & operator >> (std::istream & in, QuadraticIdeal<T> & A);

  template <class T>
  std::ostream & operator << (std::ostream & out, const QuadraticIdeal<T> & A);


  //
  //
  //
  // And now, the class
  //
  //
  // (num/dem)(aZ + (b+sqrt(D))/2 Z
  //

  template <class T>
  class QuadraticIdeal: public QuadraticIdealBase<T>
  {
  protected:
    QQ<T> q; // (num/dem) coefficient
		
  public:
    QuadraticIdeal();
    QuadraticIdeal(const QuadraticIdeal<T> & quad_i);
    ~QuadraticIdeal();

    QuadraticIdeal(const QuadraticIdealBase<T> & quad_f);
    QuadraticIdeal(const QuadraticNumber<T> & quad_n);
    QuadraticIdeal(const QQ<T> & q_in, const QuadraticIdealBase<T> & quad_f);
    QuadraticIdeal(const T & num, const T& den,
		    const QuadraticIdealBase<T> & quad_f);
    QuadraticIdeal(const T & num, const T& den,
		    const T & a_in, const T & b_in);

    void assign_one();
    bool assign_prime(const T & p);
    //		bool assign_prime(const long p);
    /*		{
      return assign_prime(to_ZZ(p));
      }
    */
    bool assign(const T & a_in, const T & b_in);
    void assign(const QuadraticIdealBase<T> & quad_f);
    bool assign(const QuadraticNumber<T> & quad_n);
    bool assign(const QQ<T> & q_in, const QuadraticIdealBase<T> & quad_f);
    bool assign(const T & num, const T & den, const QuadraticIdealBase<T> & quad_f);
    bool assign(const QQ<T> & q_in, const T & a_in, const T & b_in);
    bool assign(const T & num, const T & den, const T & a_in, const T & b_in);
    void assign(const QuadraticIdeal<T> & quad_i);

    /// Sets \p n and \p d to the numerator and denominator.
    void get_q(T& n, T& d);
    void get_q(QQ<T> & q_out);
    void set_q(const QQ<T> & q_in);
    void set_q(const T & n);
    void set_inv_q(const T & d);

    void reduce_imag();
    void reduce_imag(QuadraticNumber<T> &quad_n);
    void reduce_real();
    void reduce_real(QuadraticNumber<T> &quad_n);
    void reduce_real(vector< QuadraticNumber<T> > &quad_n);
    void reduce();

    /**
     * @brief Returns true iff the conditions in [Jac99, sec 2.1] are
     * satisfied.
     */
    bool is_standard_representation() const;

    /**
     * @brief Returns true iff the conditions in [Jac99, sec 2.1] are
     * satisfied.
     */
    bool is_normal() const;

    /**
     * @brief Returns true iff [Jac99, def 2.11] is satisfied.
     *
     * The ideal is reduced if it is of the form
     * I = 1/a (Z a + Z (b+sqrt(D))/2) with (a,b,c) reduced, which is
     * equivalent to the fact that 1 is a minimum in I.
     * 
     * This ideal must be in standard form and must be normal, or an
     * assertion will fail.
     */
    bool is_reduced() const;

    friend void multiply<T>(QuadraticIdeal<T> & C,
			    const QuadraticIdeal<T> & A, const QuadraticIdeal<T> &B);

    friend void multiply_real<T>(QuadraticIdeal<T> & C,
				 const QuadraticIdeal<T> & A,
				 const QuadraticIdeal<T> & B,
				 QuadraticNumber<T> & gamma,
				 const QuadraticNumber<T> & alpha,
				 const QuadraticNumber<T> & beta);

		
    friend void square<T>(QuadraticIdeal<T> & C,
			  const QuadraticIdeal<T> & A);


    friend void power<T> (QuadraticIdeal<T> &C,
			  const QuadraticIdeal<T> & A, const ZZ & n);

    friend QuadraticIdeal<T> inverse<T> (const QuadraticIdeal<T>
					  &A);

    friend QuadraticIdeal<T> operator - <T>(const QuadraticIdeal<T> & A);

    friend QuadraticIdeal<T> operator * <T>(const QuadraticIdeal<T> & A,
					     const QuadraticIdeal<T> & B);

    QuadraticIdeal<T> & operator *= (const QuadraticIdeal<T> & A);


    void norm(T & n, T & d);
    void norm(QQ<T> & nm);



    friend void swap<T> (QuadraticIdeal<T> & A, QuadraticIdeal<T> & B);



    bool is_one() const;
    bool is_equal(const QuadraticIdeal<T> & B) const;

    friend bool operator == <T> (const QuadraticIdeal<T> & A,
				 const QuadraticIdeal<T> & B);

    friend bool operator != <T> (const QuadraticIdeal<T> & A,
				 const QuadraticIdeal<T> & B);

    friend std::istream & operator >> <T> (std::istream & in,
					   QuadraticIdeal<T> & A);
    friend std::ostream & operator << <T> (std::ostream & out,
					   const QuadraticIdeal<T> & A);

  };

  //
  // Declare Specialized Methods
  //
  template <> QuadraticIdeal<ZZ>::QuadraticIdeal(const QuadraticNumber<ZZ> & quad_n);
  //	template <> bool QuadraticIdeal<long>::assign_prime(const long p);
  //	template <> bool QuadraticIdeal<ZZ>::assign_prime(const long p);
  template <> bool QuadraticIdeal<ZZ>::assign(const QuadraticNumber<ZZ> & quad_n);
  template <> bool QuadraticIdeal<ZZ>::is_normal() const;
  template <> bool QuadraticIdeal<ZZ>::is_standard_representation() const;
  template <> bool QuadraticIdeal<ZZ>::is_reduced() const;
  template <> bool QuadraticIdeal<long long>::assign(const QuadraticNumber<long long> & quad_n);
  template <> bool QuadraticIdeal<GF2EX>::assign(const QuadraticNumber<GF2EX> & quad_n);
  template <> QuadraticIdeal<GF2EX> inverse(const QuadraticIdeal<GF2EX> &A);


} // ANTL

// Unspecialized template definitions.
//#include "../../../src/quadratic/QuadraticIdeal_impl.hpp"

#endif // guard

