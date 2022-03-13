/**
 * @file HashEntry_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from HashEntry.h only.
 * @version $Header$
 */

namespace ANTL {

template <class T>
bool operator==(const HashEntry<T> &A, const HashEntry<T> &B) {
  if (A.a != B.a)
    return false;
  return (A.b == B.b);
}

template <>
inline bool operator==<ZZ>(const HashEntry<ZZ> &A, const HashEntry<ZZ> &B) {
  if (A.a != B.a)
    return false;
  return ((A.b % (A.a << 1)) == (B.b % (B.a << 1)));
}

template <>
inline bool operator==
    <long long>(const HashEntry<long long> &A, const HashEntry<long long> &B) {
  if (A.a != B.a)
    return false;
  return ((A.b % (A.a << 1)) == (B.b % (B.a << 1)));
}

template <>
inline bool operator==
    <long>(const HashEntry<long> &A, const HashEntry<long> &B) {
  if (A.a != B.a)
    return false;
  return ((A.b % (A.a << 1)) == (B.b % (B.a << 1)));
}

template <class T>
bool operator!=(const HashEntry<T> &A, const HashEntry<T> &B) {
  return (!(A == B));
}

//
// operator >>
//
// Task:
//      inputs a qi_pair<T> from the istream in.
//

template <class T> std::istream &operator>>(std::istream &in, HashEntry<T> &A) {
  long n = 0;
  char c;
  T ibuf[2];

  in >> c;
  if (c != '(') {
    cout << "ERROR:  HashEntry::operator>>::( expected" << endl;
    exit(1);
  }

  in >> c;
  while (c != ')' && n != 2) {
    in.putback(c);
    in >> ibuf[n];
    n++;
    in >> c;
    if (c == ',')
      in >> c;
  }

  A.a = ibuf[0];
  A.b = ibuf[1];

  return in;
}

//
// operator <<
//
// Task:
//      outputs a qi_pair<T> to the ostream out.
//

template <class T>
std::ostream &operator<<(std::ostream &out, const HashEntry<T> &A) {
  out << "(" << A.a << ", " << A.b << ")" << flush;
  return out;
}

} // namespace ANTL
