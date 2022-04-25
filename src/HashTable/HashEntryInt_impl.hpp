/**
 * @file HashEntryInt_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from HashEntryInt.h only.
 * @version $Header$
 */

namespace ANTL {

template <class T, class S>
HashEntryInt<T, S>::HashEntryInt(const T &na, const T &nb, const S &nd)
    : HashEntry<T>(na, nb) {
  d = nd;
}

//
// operator >>
//
// Task:
//      inputs a qi_pair<T> from the istream in.
//

template <class T, class S>
std::istream &operator>>(std::istream &in, HashEntryInt<T, S> &A) {
  long n = 0;
  char c;
  T ibuf[2];
  S dist;

  in >> c;
  if (c != '(') {
    cout << "ERROR:  HashEntryInt::operator>>::( expected" << endl;
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
std::ostream &operator<<(std::ostream &out, const HashEntryInt<T, S> &A) {
  out << "(" << A.a << ", " << A.b << ", " << A.d << ")" << endl;
  return out;
}

} // namespace ANTL
