/**
 * @file qo_hash_entry_real.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from HashEntryVec.h only.
 * @version $Header$
 */

namespace ANTL {

template <class T>
HashEntryVec<T>::HashEntryVec(const T &na, const T &nb, const vec_ZZ &nd)
    : HashEntry<T>(na, nb) {
  d = nd;
}

//
// operator >>
//
// Task:
//      inputs a qi_pair<T> from the istream in.
//

template <class T>
std::istream &operator>>(std::istream &in, HashEntryVec<T> &A) {
  long n = 0;
  char c;
  T ibuf[2];
  vec_ZZ dist;

  in >> c;
  if (c != '(') {
    cout << "ERROR:  HashEntryVec::operator>>::( expected" << endl;
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

template <class T>
std::ostream &operator<<(std::ostream &out, const HashEntryVec<T> &A) {
  out << "(" << A.a << ", " << A.b << ", " << A.d << ")" << endl;
  return out;
}

} // namespace ANTL
