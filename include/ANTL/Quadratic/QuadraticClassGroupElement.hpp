/**
 * @file QuadraticClassGroupElement.hpp
 * @author Michael Jacobson (MJ)
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_QI_CLASS_H
#define ANTL_QUADRATIC_QI_CLASS_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

namespace ANTL
{

  template < class T > class QuadraticClassGroupElement;

  template < class T >
  void multiply_imag(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A,
		     const QuadraticClassGroupElement<T> & B);

  template < class T >
  void multiply_real(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A,
		     const QuadraticClassGroupElement<T> & B);

  template < class T >
  void multiply(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A,
		const QuadraticClassGroupElement<T> & B);

  template < class T >
  void square_imag(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A);

  template < class T >
  void square_real(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A);

  template < class T >
  void square(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A);

  template < class T >
  void power_imag(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A, const ZZ & n);

  template < class T >
  void power_real(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A, const ZZ & n);

  template < class T >
  void power(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A, const ZZ & n);

  template < class T >
  void nupower_imag(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A, const ZZ & n);

  template < class T >
  void nupower_real(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A, const ZZ & n);

  template < class T >
  void nupower(QuadraticClassGroupElement<T> & C, const QuadraticClassGroupElement<T> & A, const ZZ & n);

  template < class T >
  QuadraticClassGroupElement<T> operator * (const QuadraticClassGroupElement<T> & A, const QuadraticClassGroupElement<T> & B);

  template < class T >
  void swap(QuadraticClassGroupElement<T> & A, QuadraticClassGroupElement<T> & B);

  template < class T >
  std::istream & operator >> (std::istream & in, QuadraticClassGroupElement<T> & A);

  //
  // Class: QuadraticClassGroupElement<T>
  //
  // This class represents an element in the class group of a quadratic
  //    order, i.e., a reduced, primitive, invertible ideal
  //    aZ + (b+sqrt(Delta))/2 Z  where Delta is the discriminant of the
  //    quadratic order.  This class is derived from QuadraticIdealBase and shares all
  //    its properties, but each instance is forced to be reduced.
  //

  template < class  T >
  class QuadraticClassGroupElement : public QuadraticIdealBase<T>
  {
  public:

    //
    // constructors and destructor
    //

    QuadraticClassGroupElement();
    QuadraticClassGroupElement(const QuadraticIdealBase<T> & A);
    QuadraticClassGroupElement(const HashEntry<T> & A);

    template < class S >
    QuadraticClassGroupElement(const HashEntryInt<T,S> & A);

    QuadraticClassGroupElement(const HashEntryVec<T> & A);
    QuadraticClassGroupElement(const QuadraticClassGroupElement<T> & A);

    ~QuadraticClassGroupElement();


    /*
    //
    // initialization
    //

    static void set_current_order(const T & newf, const T & newh);
    static void set_current_order(const T & newDelta);
    */


    //
    // assignment
    //

    bool assign_prime(const T & p);
    bool assign(const T & a2, const T & b2);
    void assign(const HashEntry<T> & B);

    template < class S >
    void assign(const HashEntryInt<T,S> & B);

    void assign(const HashEntryVec<T> & B);
    void assign(const QuadraticClassGroupElement<T> & B);
    QuadraticClassGroupElement<T> & operator = (const QuadraticClassGroupElement<T> & A);


    //
    // arithmetic operations
    //

    friend void multiply_imag <T> (QuadraticClassGroupElement<T> & C,
                                   const QuadraticClassGroupElement<T> & A,
                                   const QuadraticClassGroupElement<T> & B);

    friend void multiply_real <T> (QuadraticClassGroupElement<T> & C,
                                   const QuadraticClassGroupElement<T> & A,
                                   const QuadraticClassGroupElement<T> & B);

    friend void multiply <T> (QuadraticClassGroupElement<T> & C,
			      const QuadraticClassGroupElement<T> & A,
			      const QuadraticClassGroupElement<T> & B);

    friend void square_imag <T> (QuadraticClassGroupElement<T> & C,
                                 const QuadraticClassGroupElement<T> & A);
    friend void square_real <T> (QuadraticClassGroupElement<T> & C,
                                 const QuadraticClassGroupElement<T> & A);
    friend void square <T> (QuadraticClassGroupElement<T> & C,
			    const QuadraticClassGroupElement<T> & A);

    friend void power_imag <T> (QuadraticClassGroupElement<T> & C,
				const QuadraticClassGroupElement<T> & A,
				const ZZ & n);
    friend void power_real <T> (QuadraticClassGroupElement<T> & C,
				const QuadraticClassGroupElement<T> & A,
				const ZZ & n);
    friend void power <T> (QuadraticClassGroupElement<T> & C,
			   const QuadraticClassGroupElement<T> & A,
			   const ZZ & n);

    friend void nupower_imag <T> (QuadraticClassGroupElement<T> & C,
				  const QuadraticClassGroupElement<T> & A,
				  const ZZ & n);
    friend void nupower_real <T> (QuadraticClassGroupElement<T> & C,
				  const QuadraticClassGroupElement<T> & A,
				  const ZZ & n);
    friend void nupower <T> (QuadraticClassGroupElement<T> & C,
			     const QuadraticClassGroupElement<T> & A,
			     const ZZ & n);

    friend QuadraticClassGroupElement<T> operator * <T> (const QuadraticClassGroupElement<T> & A,
				       const QuadraticClassGroupElement<T> & B);
    QuadraticClassGroupElement<T> & operator *= (const QuadraticClassGroupElement<T> & A);



    //
    // basic functions
    //

    friend void swap <T> (QuadraticClassGroupElement<T> & A, QuadraticClassGroupElement<T> & B);
    friend ZZ hash_value (const QuadraticClassGroupElement<T> &A)
    { return hash_value_base(A.a); }




    //
    // input/output
    //

    friend std::istream & operator >> <T> (std::istream & in, QuadraticClassGroupElement<T> & A);
  };

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Quadratic/QuadraticClassGroupElement_impl.hpp"

#endif // guard

