/**
 * @file quadratic_ideal.hpp
 * @version $Header: /scrinium/ANTL/ANTL/include/ANTL/quadratic/quadratic_ideal.hpp,v 1.1 2013/02/26 06:29:10 jacobs Exp $
 *
 * Basically a quadratic form, but with additional q, so that the ideal is
 * I = (num/dem)(aZ + ((b+Sqrt(Delta))/2) Z)
 */

#ifndef ANTL_QUADRATIC_QUADRATIC_IDEAL_H
#define ANTL_QUADRATIC_QUADRATIC_IDEAL_H

#include <ANTL/Quadratic/QuadraticForm.hpp>
#include <ANTL/Arithmetic/QQ.hpp>

namespace ANTL
{

  template < class T> class quadratic_ideal;

  template <class T>
  void multiply(quadratic_ideal<T> & C,
		const quadratic_ideal<T> & A, const quadratic_ideal<T> &B);

  template < class T >
  void multiply_real(quadratic_ideal<T> & C, const quadratic_ideal<T> & A,
		     const quadratic_ideal<T> & B, quadratic_number<T> & gamma,
		     const quadratic_number<T> & alpha, const quadratic_number<T> & beta);

  //similarly, needed for squaring, powers
  template <class T>
  void square(quadratic_ideal<T> & C, const quadratic_ideal<T> & A);

  template <class T>
  void power (quadratic_ideal<T> &C, const quadratic_ideal<T> & A, const ZZ & n);


  template <class T>
  quadratic_ideal<T> inverse(const quadratic_ideal<T> & A);

  template <class T>
  quadratic_ideal<T> operator - (const quadratic_ideal<T> & A);

  template <class T>
  quadratic_ideal<T> operator * (const quadratic_ideal<T> & A,
				 const quadratic_ideal<T> & B);

  template <class T>
  void swap (quadratic_ideal<T> & A, quadratic_ideal<T> & B);

  template <class T>
  bool operator == (const quadratic_ideal<T> & A, const quadratic_ideal<T> & B);

  template <class T>
  bool operator != (const quadratic_ideal<T> & A, const quadratic_ideal<T> & B);

  template <class T>
  std::istream & operator >> (std::istream & in, quadratic_ideal<T> & A);

  template <class T>
  std::ostream & operator << (std::ostream & out, const quadratic_ideal<T> & A);


  //
  //
  //
  // And now, the class
  //
  //
  // (num/dem)(aZ + (b+sqrt(D))/2 Z
  //

  template <class T>
  class quadratic_ideal: public quadratic_form<T>
  {
  protected:
    QQ<T> q; // (num/dem) coefficient
		
  public:
    quadratic_ideal();
    quadratic_ideal(const quadratic_ideal<T> & quad_i);
    ~quadratic_ideal();

    quadratic_ideal(const quadratic_form<T> & quad_f);
    quadratic_ideal(const quadratic_number<T> & quad_n);
    quadratic_ideal(const QQ<T> & q_in, const quadratic_form<T> & quad_f);
    quadratic_ideal(const T & num, const T& den,
		    const quadratic_form<T> & quad_f);
    quadratic_ideal(const T & num, const T& den,
		    const T & a_in, const T & b_in);

    void assign_one();
    bool assign_prime(const T & p);
    //		bool assign_prime(const long p);
    /*		{
      return assign_prime(to_ZZ(p));
      }
    */
    bool assign(const T & a_in, const T & b_in);
    void assign(const quadratic_form<T> & quad_f);
    bool assign(const quadratic_number<T> & quad_n);
    bool assign(const QQ<T> & q_in, const quadratic_form<T> & quad_f);
    bool assign(const T & num, const T & den, const quadratic_form<T> & quad_f);
    bool assign(const QQ<T> & q_in, const T & a_in, const T & b_in);
    bool assign(const T & num, const T & den, const T & a_in, const T & b_in);
    void assign(const quadratic_ideal<T> & quad_i);

    /// Sets \p n and \p d to the numerator and denominator.
    void get_q(T& n, T& d);
    void get_q(QQ<T> & q_out);
    void set_q(const QQ<T> & q_in);
    void set_q(const T & n);
    void set_inv_q(const T & d);

    void reduce_imag();
    void reduce_imag(quadratic_number<T> &quad_n);
    void reduce_real();
    void reduce_real(quadratic_number<T> &quad_n);
    void reduce_real(vector< quadratic_number<T> > &quad_n);
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

    friend void multiply<T>(quadratic_ideal<T> & C,
			    const quadratic_ideal<T> & A, const quadratic_ideal<T> &B);

    friend void multiply_real<T>(quadratic_ideal<T> & C,
				 const quadratic_ideal<T> & A,
				 const quadratic_ideal<T> & B,
				 quadratic_number<T> & gamma,
				 const quadratic_number<T> & alpha,
				 const quadratic_number<T> & beta);

		
    friend void square<T>(quadratic_ideal<T> & C,
			  const quadratic_ideal<T> & A);


    friend void power<T> (quadratic_ideal<T> &C,
			  const quadratic_ideal<T> & A, const ZZ & n);

    friend quadratic_ideal<T> inverse<T> (const quadratic_ideal<T>
					  &A);

    friend quadratic_ideal<T> operator - <T>(const quadratic_ideal<T> & A);

    friend quadratic_ideal<T> operator * <T>(const quadratic_ideal<T> & A,
					     const quadratic_ideal<T> & B);

    quadratic_ideal<T> & operator *= (const quadratic_ideal<T> & A);


    void norm(T & n, T & d);
    void norm(QQ<T> & nm);



    friend void swap<T> (quadratic_ideal<T> & A, quadratic_ideal<T> & B);



    bool is_one() const;
    bool is_equal(const quadratic_ideal<T> & B) const;

    friend bool operator == <T> (const quadratic_ideal<T> & A,
				 const quadratic_ideal<T> & B);

    friend bool operator != <T> (const quadratic_ideal<T> & A,
				 const quadratic_ideal<T> & B);

    friend std::istream & operator >> <T> (std::istream & in,
					   quadratic_ideal<T> & A);
    friend std::ostream & operator << <T> (std::ostream & out,
					   const quadratic_ideal<T> & A);

  };

  //
  // Declare Specialized Methods
  //
  template <> quadratic_ideal<ZZ>::quadratic_ideal(const quadratic_number<ZZ> & quad_n);
  //	template <> bool quadratic_ideal<long>::assign_prime(const long p);
  //	template <> bool quadratic_ideal<ZZ>::assign_prime(const long p);
  template <> bool quadratic_ideal<ZZ>::assign(const quadratic_number<ZZ> & quad_n);
  template <> bool quadratic_ideal<ZZ>::is_normal() const;
  template <> bool quadratic_ideal<ZZ>::is_standard_representation() const;
  template <> bool quadratic_ideal<ZZ>::is_reduced() const;
  template <> bool quadratic_ideal<long long>::assign(const quadratic_number<long long> & quad_n);
  template <> bool quadratic_ideal<GF2EX>::assign(const quadratic_number<GF2EX> & quad_n);
  template <> quadratic_ideal<GF2EX> inverse(const quadratic_ideal<GF2EX> &A);


} // ANTL

// Unspecialized template definitions.
//#include "../../../src/quadratic/quadratic_ideal_impl.hpp"

#endif // guard

