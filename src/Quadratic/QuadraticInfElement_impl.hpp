namespace ANTL {

template <class T>
QuadraticInfElement<T>::QuadraticInfElement(QuadraticOrder<T> & quad_o) : qib(quad_o) {
  qib.assign_one();
  Delta = quad_o.getDiscriminant();
  FloorRootDelta = FloorToZZ(sqrt(to_RR(Delta)));
  distance = 1;
}

// destructor
template <class T>
QuadraticInfElement<T>::~QuadraticInfElement () {}

template <class T>
QuadraticIdealBase<T> QuadraticInfElement<T>::get_qib(){
  return qib;
}

template <class T>
RR QuadraticInfElement<T>::get_distance(){
  return distance;
}


template <class T>
void QuadraticInfElement<T>::baby_step() {


}

template <class T>
void QuadraticInfElement<T>::giant_step(const QuadraticIdealBase<T> & quad_ib) {

}
}
