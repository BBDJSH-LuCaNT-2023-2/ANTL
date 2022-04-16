/**
 * @file HashEntryReal_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from HashEntryReal.h only.
 * @version $Header$
 */

namespace ANTL {

template <class T, class S>
HashEntryReal<T, S>::HashEntryReal(const T &na, const T &nb, const S &nd)
    : HashEntry<T>(na, nb) {
  d = nd;
}

// template <class T>
// HashEntryReal<T>::HashEntryReal(const T &na, const T &nb,
//                                 const qo_distance<T> &nd)
//     : HashEntry<T>(na, nb) {
//   d = nd;
// }

//
// operator >>
//
// Task:
//      inputs a qi_pair<T> from the istream in.
//

template <class T, class S>
std::istream &operator>>(std::istream &in, HashEntryReal<T, S> &A) {
  long n = 0;
  char c;
  T ibuf[2];

  RR dist;
  //   //qo_distance<T> dist;

  in >> c;
  if (c != '(') {
    cout << "ERROR:  HashEntryReal::operator>>::( expected" << endl;
    exit(1);
  }

  in >> c;
  while (c != ')' && n != 3) {
    in.putback(c);
    if (n == 2)
      in >> dist;
    else
      in >> ibuf[n];
    n++;
    in >> c;
    if (c == ',')
      in >> c;
  }

  A.a = ibuf[0];
  A.b = ibuf[1];
  A.d = dist;

  return in;
}

//
// operator <<
//
// Task:
//      outputs a qi_pair<T> to the ostream out.
//

template <class T, class S>
std::ostream &operator<<(std::ostream &out, const HashEntryReal<T, S> &A) {
  out << "(" << A.a << ", " << A.b << ", " << A.d << ")" << flush;
  return out;
}

} // namespace ANTL
