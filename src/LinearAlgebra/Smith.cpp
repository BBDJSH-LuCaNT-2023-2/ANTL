/**
 * @file smith.cpp
 * @author Kris Luttmer
 * @version $Header$
 */

#include <ANTL/LinearAlgebra/Smith.hpp>

namespace ANTL {
mat_ZZ SmithNormalForm(mat_ZZ A, mat_ZZ &U, mat_ZZ &V) {

  ZZ detA; // the determinant of the matrix A
  long n;  // the dimension of the matrix A
  long i, c;
  long j, k, l;        // looping variables
  ZZ b;                // a varible used to test when we are done
  ZZ dd, u, v;         // values used when computing the extended gcd
  vec_ZZ temp1, temp2; // a temporary vector to hold columns and rows
  mat_ZZ D;            // the diagonal entries of the smith normal form matrix
  ZZ r;                // used as a remainder after division
  bool done = false;   // to signify completion or not

  // check to be sure the matrix is square
  if (A.NumCols() != A.NumRows()) {
    cerr << "Cannot convert non-square matrix to Smith Normal Form" << endl;
    exit(0);
  }

  // check to be sure the matrix is nonempty
  if (A.NumCols() == 0) {
    cerr << "Cannot convert empty matrix to Smith Normal Form" << endl;
    exit(0);
  }

  if (A.NumRows() == 1) {
    D = A;
    ident(U, 1);
    ident(V, 1);
    return D;
  }

  // compute the determinant of the matrix A
  //   determinant (detA, A);
  set(detA);
  for (i = 1; i <= A.NumRows(); ++i)
    detA *= A(i, i);

  detA = abs(detA);
  // check to be sure the matrix is non-singular
  if (IsZero(detA)) {
    cerr << "Cannot convert singular matrix to Smith Normal Form" << endl;
    exit(0);
  }

  n = A.NumCols();
  i = n;
  temp1.SetLength(n);
  temp2.SetLength(n);
  D.SetDims(n, n);
  ident(U, n);
  ident(V, n);

  while (done == false) {
    do {

      // initialize j and c for row reduction
      c = 0;
      j = i;
      while (j != 1) {
        j = j - 1;

        // notice that A is being indexed from 1 not from 0
        if (A(i, j) == 0)
          continue;

        // start step 4 of the algorithm
        if ((A(i, j) % A(i, i)) == 0) {
          set(u);
          clear(v);
          dd = A(i, i);
        } else
          XGCD(dd, u, v, A(i, i), A(i, j));

        // use r and b as temporary variables to avoid
        // excessive division in the for loop.
        r = A(i, i) / dd;
        b = A(i, j) / dd;

        for (k = 1; k <= n; k++) {
          temp1(k) = (u * A(k, i) + v * A(k, j));
          temp2(k) = (u * V(k, i) + v * V(k, j));
          A(k, j) = (r * A(k, j) - b * A(k, i));
          V(k, j) = (r * V(k, j) - b * V(k, i));
          A(k, i) = temp1(k);
          V(k, i) = temp2(k);
        }
      }
#ifdef DEBUG
      cout << "Done Transformation" << endl;
      cout << "dd = " << dd << ", u = " << u << ", v = " << v << endl;
      cout << A << endl << V << endl;
#endif
      // initialize j and c for column reduction
      j = i;
      while (j != 1) {
        j = j - 1;

        // notice that A is being indexed from 1 not from 0
        if (A(j, i) == 0)
          continue;

        // start step 7 of the algorithm
        if ((A(j, i) % A(i, i)) == 0) {
          set(u);
          clear(v);
          dd = A(i, i);
        } else
          XGCD(dd, u, v, A(i, i), A(j, i));

        temp1 = u * A(i) + v * A(j);
        temp2 = u * U(i) + v * U(j);
        U(j) = (A(i, i) / dd) * U(j) - (A(j, i) / dd) * U(i);
        A(j) = (A(i, i) / dd) * A(j) - (A(j, i) / dd) * A(i);
        A(i) = temp1;
        U(i) = temp2;
        // make sure all the numbers are modulo the determinant
        /*  for(k = 1; k <= n; k++)
        {A(j,k) %= detA; A(i,k) %= detA;} */
        // U(i,k) %= detA; U(j,k) %= detA;}

#ifdef DEBUG
        cout << "Done Transform 2" << endl;
        cout << "dd = " << dd << ", u = " << u << ", v = " << v << endl;
        cout << A << endl << V << endl;
#endif
        c++;
      }
    } while (c > 0);

    // now do step 9
    b = A(i, i);
    done = true;

    for (k = 1; k < i && done == true; k++) {
      for (l = 1; l < i && done == true; l++) {
        rem(r, A(k, l), b);
        if (!IsZero(r)) {
          A(i) = A(i) + A(k);
          U(i) = U(i) + U(k);
          done = false;
        }
      }
    }

#ifdef DEBUG
    cout << "Done Step 9" << endl;
    cout << A << endl << U << endl;
#endif

    if (done == true) {
      done = false;
      // do step 10 here
      D(i, i) = GCD(A(i, i), detA);
      detA = detA / D(i, i);

      if (i == 2)
        done = true, D(1, 1) = GCD(A(1, 1), detA);
      i = i - 1;
    }
  }

  for (i = 1; i <= n; i++) {
    if (A(i, i) < 0) {
      U(i) *= -1;
      A(i) *= -1;
    }
  }
  return D;
}
}// ANTL
