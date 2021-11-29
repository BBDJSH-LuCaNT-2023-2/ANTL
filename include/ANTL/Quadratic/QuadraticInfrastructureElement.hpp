/**
 * @file quadratic_ideal_pair.hpp
 * @version $Header: /scrinium/ANTL/ANTL/include/ANTL/quadratic/quadratic_ideal_pair.hpp,v 1.1 2013/02/26 06:29:10 jacobs Exp $
 *
 * A quadratic ideal (not necessarily reduced) and relative generator of type
 * quadratic number
 */

#ifndef ANTL_QUADRATIC_IDEAL_PAIR_H
#define ANTL_QUADRATIC_IDEAL_PAIR_H

#include <ANTL/quadratic/quadratic_form.hpp>
#include <ANTL/quadratic/qi_class.hpp>
#include <ANTL/quadratic/quadratic_number.hpp>

namespace ANTL
{

  template < class T> class quadratic_ideal_pair;

  template <class T>
  void multiply(quadratic_ideal_pair<T> & C,
		const quadratic_ideal_pair<T> & A, const quadratic_ideal_pair<T> &B);

  template <class T>
  void square(quadratic_ideal_pair<T> & C, const quadratic_ideal_pair<T> & A);

  template <class T>
  void power (quadratic_ideal_pair<T> &C, const quadratic_ideal_pair<T> & A, const ZZ & n);

  template <class T>
  void power_reduce (quadratic_ideal_pair<T> &C, const quadratic_ideal_pair<T> & A, const ZZ & n);

  template <class T>
  quadratic_ideal_pair<T> inverse(const quadratic_ideal_pair<T> & A);

  template <class T>
  quadratic_ideal_pair<T> operator - (const quadratic_ideal_pair<T> & A);

  template <class T>
  quadratic_ideal_pair<T> operator * (const quadratic_ideal_pair<T> & A,
				      const quadratic_ideal_pair<T> & B);

  template <class T>
  void swap (quadratic_ideal_pair<T> & A, quadratic_ideal_pair<T> & B);

  template <class T>
  bool operator == (const quadratic_ideal_pair<T> & A, const quadratic_ideal_pair<T> & B);

  template <class T>
  bool operator != (const quadratic_ideal_pair<T> & A, const quadratic_ideal_pair<T> & B);

  template <class T>
  std::istream & operator >> (std::istream & in, quadratic_ideal_pair<T> & A);

  template <class T>
  std::ostream & operator << (std::ostream & out, const quadratic_ideal_pair<T> & A);


  //
  // And now, the class
  //
  //
  // (alpha)(aZ + (b+sqrt(D))/2 Z
  //

  template <class T>
  class quadratic_ideal_pair: public quadratic_form<T>
  {
  protected:
    quadratic_number<T> alpha;

  public:
    quadratic_ideal_pair();
    ~quadratic_ideal_pair();

    //
    // assignments
    //

    void assign_one();
    bool assign_prime(const T & p);
    bool assign(const T & a_in, const T & b_in);
    void assign(const quadratic_form<T> & q_f);
    void assign(const qi_class<T> & q_i);
    bool assign(const quadratic_number<T> & q_n);
    void assign(const quadratic_number<T> & q_n,
		const quadratic_form<T> & q_f);
    void assign(const quadratic_number<T> & q_n,
		const qi_class<T> & q_i);
    void assign(const quadratic_ideal_pair<T> & q_i);


    //
    // access
    //

    void get_alpha(quadratic_number<T> & a)  { a.assign(alpha); }
    void set_alpha(quadratic_number<T> & a) {alpha = a;}


    //
    // reduction
    //

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


    //
    // arithmetic
    //

    void rho();
    void inverse_rho();

    friend void multiply<T>(quadratic_ideal_pair<T> & C,
			    const quadratic_ideal_pair<T> & A,
			    const quadratic_ideal_pair<T> &B);
		
    friend void square<T>(quadratic_ideal_pair<T> & C,
			  const quadratic_ideal_pair<T> & A);

    friend void power<T> (quadratic_ideal_pair<T> &C,
			  const quadratic_ideal_pair<T> & A,
			  const ZZ & n);

    friend void power_reduce<T> (quadratic_ideal_pair<T> &C,
				 const quadratic_ideal_pair<T> & A,
				 const ZZ & n);

    friend quadratic_ideal_pair<T> inverse<T> (const quadratic_ideal_pair<T> &A);

    friend quadratic_ideal_pair<T> operator - <T>(const quadratic_ideal_pair<T> & A);

    friend quadratic_ideal_pair<T> operator * <T>(const quadratic_ideal_pair<T> & A,
						  const quadratic_ideal_pair<T> & B);

    quadratic_ideal_pair<T> & operator *= (const quadratic_ideal_pair<T> & A);


    T norm();



    friend void swap<T> (quadratic_ideal_pair<T> & A,
			 quadratic_ideal_pair<T> & B);


    //
    // comparisons
    //

    bool is_one() const;
    bool is_equal(const quadratic_ideal_pair<T> & B) const;
    bool is_equal(const quadratic_number<T> & qn) const;

    friend bool operator == <T> (const quadratic_ideal_pair<T> & A,
				 const quadratic_ideal_pair<T> & B);

    friend bool operator != <T> (const quadratic_ideal_pair<T> & A,
				 const quadratic_ideal_pair<T> & B);

    friend std::istream & operator >> <T> (std::istream & in,
					   quadratic_ideal_pair<T> & A);
    friend std::ostream & operator << <T> (std::ostream & out,
					   const quadratic_ideal_pair<T> & A);

  };

  //
  // Declare Specialized Methods
  //

  template <> bool quadratic_ideal_pair<ZZ>::assign(const quadratic_number<ZZ> & q_n);
  template <> bool quadratic_ideal_pair<ZZ>::is_normal() const;
  template <> bool quadratic_ideal_pair<ZZ>::is_standard_representation() const;
  template <> bool quadratic_ideal_pair<ZZ>::is_reduced() const;
  template <> bool quadratic_ideal_pair<long long>::assign(const quadratic_number<long long> & q_n);
  template <> bool quadratic_ideal_pair<GF2EX>::assign(const quadratic_number<GF2EX> & q_n);
  template <> quadratic_ideal_pair<GF2EX> inverse(const quadratic_ideal_pair<GF2EX> &A);


} // ANTL

// Unspecialized template definitions.
#include "../../../src/quadratic/quadratic_ideal_pair_impl.hpp"

#endif // guard
