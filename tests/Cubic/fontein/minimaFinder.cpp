#include "../cubic_header.hpp"
#include <vector>

int main(){

  ZZ real_ibcf[4];
  real_ibcf[3] = 2;
  real_ibcf[2] = 1;
  real_ibcf[1] = -6;
  real_ibcf[0] = -2;

  polynomial<ZZ> const real_poly{{real_ibcf[0],real_ibcf[1],real_ibcf[2],real_ibcf[3] }};
  CubicOrder<ZZ, RR> * ro_point; ro_point = CubicOrder<ZZ, RR>::make_order(real_poly);
  CubicOrder<ZZ, RR> * Odessa = ro_point;

  cout << "rho1:  " << Odessa->get_rho1() << "     rho2:  " << Odessa->get_rho2() << endl;
  cout << " root1 " << Odessa->get_root1() << " root2 " << Odessa->get_root2() << " root3 " << Odessa->get_root3() << std::endl;
  cout << "order disc:  "<< Odessa->get_discriminant() << endl;
  cout << "----------------------------------------" << endl;

  int xctr = 0;

  RR logholder, rtemp;
  RR xbound, zbound, ybound;
  xbound = RR(15);
  zbound = RR(15);
  ybound = RR(15);

  QQ<ZZ> normie;
  logholder = RR(0);
  std::vector< std::vector<RR> > logarray;
  std::vector< std::vector<CubicIdeal<ZZ,RR>> > ideal_list;

  std::vector<CubicIdeal<ZZ,RR>> ideal_list_l2;
  std::vector<RR> logarray_l2;

  CubicIdeal<ZZ,RR> starting_ideal(Odessa);
  ideal_list.push_back(ideal_list_l2);
  ideal_list[0].push_back(starting_ideal);
  CubicElement<ZZ,RR> minimum_holder(Odessa);
  logarray.push_back(logarray_l2);
  logarray[0].push_back(RR(0));

  // Compute along the X axis.
  do{
    xctr++;
    starting_ideal.adjacent_ideal(starting_ideal, minimum_holder, 'X');
    ideal_list.push_back(ideal_list_l2);
    ideal_list[xctr].push_back(starting_ideal);

    cout << "Adj element:" << minimum_holder.toString() << endl;
    minimum_holder.norm(normie);
    cout << "Adj element norm " << normie << endl;
    minimum_holder.get_real_value(rtemp);
    abs(rtemp,rtemp);
    log(rtemp,rtemp);
    add(logholder, logholder, rtemp);
    cout << "Minima's log value: " <<logholder << endl;
    logarray.push_back(logarray_l2);
    logarray[xctr].push_back(logholder);
  }while(logholder < xbound);

  for (int i =0; i < ideal_list.size(); i++){
    std::cout << ideal_list[i][0].get_coeff(0,0) << " " << ideal_list[i][0].get_coeff(1,0) << " " << ideal_list[i][0].get_coeff(2,0) << " " << std::endl;
  }


    std::cout << "Fundamental Units: " << std::endl;
  std::cout << "[" << Odessa->get_fundamental_unit(0)->get_u() << " " << Odessa->get_fundamental_unit(0)->get_x() << " " << Odessa->get_fundamental_unit(0)->get_y() << "]" <<std::endl;
  std::cout << "[" << Odessa->get_fundamental_unit(1)->get_u() << " " << Odessa->get_fundamental_unit(1)->get_x() << " " << Odessa->get_fundamental_unit(1)->get_y() << "]" <<std::endl;
  std::cout <<  "Regulator: " << Odessa->get_regulator() << std::endl;

  return 0;
}
