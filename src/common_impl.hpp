/**
 * @file utilities_impl.hpp
 * @author Michael Jacobson
 * @remark generic implementations of utility routines.
 */


//
// eval_poly - evaluate polynomial at q
//


template < class T >
ZZ eval_poly(const T & p, const ZZ & q)
{
  ZZ hval;
  static vec_ZZ qvec;
  static ZZ pp;
  long i,j,da;

  // make sure static array of powers of q has enough elements!
  if (qvec.length() == 0 || (qvec.length() > 1 && qvec[1] != q)) {
    qvec.SetLength(1);
    qvec[0] = 1;
  }

  if (qvec.length() <= deg(p)) {
    j = qvec.length();
    qvec.SetLength(deg(p)+1);

    for (i=j; i<=deg(p); ++i)
      qvec[i] = qvec[i-1]*q;
  }

  clear(hval);
  da = deg (p);
  //cout << "ZZ qvec = " << qvec << endl;
  for (i=0; i<=da; ++i)
    hval += rep(coeff(p,i)) * qvec[i];

  return hval;
}



//
// Quadratic residuosity functions
//

/*
 * Function: Jacobi_base
 * Purpose: tests to see if a is a square mod p
 */

template < class T, class U >
  inline long
  Jacobi_base (const T & a, const T & n, T & P1, T & P2, U & c, U & q)
{
  long jac = 1;

  P1 = a;
  P2 = n;

  while (deg (P1) > 0)
    {
      c = rep (LeadCoeff (P1));
      if (deg (P2) & 1)
	jac *= ::Jacobi(c, q);

      if (jac == 0)
	return jac;

      MakeMonic (P1);

      if ((deg (P1) & 1) && (deg (P2) & 1))
	jac *= ::Jacobi(q - 1, q);

      swap (P1, P2);
      P1 %= P2;
    }

  if (deg (P2) & 1)
    {
      c = rep (LeadCoeff (P1));
      jac *= ::Jacobi(c, q);
    }

  return jac;
}

//
// modular square root functions
//

template < class T >
long ressol(T & x, const T & a, const T & p)
{
  T aa,g,h,temp,temp2;
  long jac,e,s;
  long i;
  ZZ q,t;

  jac = Jacobi(a,p);

  if (jac <= 0) {
    clear(x);
    return jac;
  }

  aa = a % p;

  // compute order mod Q
  get_qdp(q,p);

  if ((q % 4) == 3)
    PowerMod(x,aa,(q+1)/4,p);
  else if ((q % 8) == 5) {
    PowerMod(temp,aa,(q-1)/4,p);
    if (IsOne(temp))
      PowerMod(x,aa,(q+3)/8,p);
    else {
      PowerMod(x,aa*4,(q-5)/8,p);
      MulMod(x,x*2,aa,p);
    }
  }
  else {
    --q;
    t = q;
    s = MakeOdd(t);

    // find quadratic non-residue
    long trys = 0;
    long pdeg = deg(p)+1;
    do {
      ++trys;
      if (trys == 10) {
        ++pdeg;
        trys = 0;
      }

      random(g,pdeg);
      g %= p;

    } while (Jacobi(g,p) != -1);

    e = 0;
    for (i=2; i<=s; ++i) {
      PowerMod(temp,InvMod(g,p),e,p);
      MulMod(temp2,temp,aa,p);
      PowerMod(temp,temp2,q >> i,p);
      if (!IsOne(temp))
        e += (1 << (i-1));
    }

    // h = a g^-e mod p
    PowerMod(temp,InvMod(g,p),e,p);
    MulMod(h,temp,aa,p);

    // x = g^(e/2) h^(t+1)/2 mod p
    PowerMod(temp2,h,(t+1) >> 1,p);
    PowerMod(temp,g,e >> 1,p);
    MulMod(x,temp2,temp,p);

    // make sure same solution is returned with every call
    T x2;
    ZZ xval, x2val;

    x2 = -x;
    x2 %= p;

    xval = eval_poly (x, ANTL::CARDINALITY<T>());
    x2val = eval_poly (x2, ANTL::CARDINALITY<T>());

    if (x2val < xval)
      x = x2;
  }

  return jac;
}


///
/// Functions used in the Cubic library
///
template<typename PType>
bool is_close(const std::vector<PType> & vec1, const std::vector<PType> & vec2, const PType & max_dist){
    if(vec1.size() != vec2.size()){
      throw "vectors are not the same size";
    }
    PType abs_dif;

    for (int i = 0; i < vec1.size(); ++i ){
      sub(abs_dif, vec1[i], vec2[i]);
      abs(abs_dif, abs_dif);
      if(abs_dif > max_dist){
        return false;
      }
    }
    return true;
};

template<typename PType>
void compute_initial_s(const std::vector<PType> & alpha, const int k_bound){
    int kprime;
    PType inf_norm = infinity_norm(alpha);
    std::vector<PType> s_term(alpha.size());
};


template<typename PType>
void log_to_valuation(std::vector<PType> &valuation_vec, const std::vector<PType> & log_vec, const int r1){
  PType tempvar, new_coord = 0;
  valuation_vec.clear();

  for(int i = 0; i < log_vec.size(); ++i ){
    if(i < r1){
      sub(new_coord, new_coord, log_vec[i]); //      new_coord -= log_vec[i];

      valuation_vec.push_back(exp(log_vec[i]));
    }else{
      mul(tempvar, 2, log_vec[i]);
      sub(new_coord, new_coord, tempvar);      //new_coord -= 2*log_vec[i];
    }
  }
  if(log_vec.size()+1 > r1){
    div(new_coord, new_coord, 2);             //new_coord = new_coord/2;

  }
  valuation_vec.push_back(exp(new_coord));

};
