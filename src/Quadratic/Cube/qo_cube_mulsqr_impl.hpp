//For some reason, the compiler doesn't like us using sqr and mul in this class   -RL

template < class T >
void
qo_cube_mulsqr<T>::cube (QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A){
    QuadraticIdealBase<T> CC(*A.get_QO());
    sqr(CC,A);
    mul(C,CC,A);
}
