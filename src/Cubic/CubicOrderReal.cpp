#ifndef ANTL_CUBIC_ORDER_REAL_CPP
#define ANTL_CUBIC_ORDER_REAL_CPP

#include "../../include/ANTL/Cubic/CubicOrderReal.hpp"



template <typename Type, typename PType>
void CubicOrderReal<Type, PType> :: compute_conjugate_elements(){

  mul(this->conjugate_bases[0][0],to<PType>((this->defining_IBCF)[3]), this->root_list[1]);      // a*delta'

  this->conjugate_bases[0][1] = this->conjugate_bases[0][0];                                          // a*delta'
  add(this->conjugate_bases[0][1],this->conjugate_bases[0][1],to<PType>((this->defining_IBCF)[2]) );  // a*delta' + b
  mul(this->conjugate_bases[0][1], this->conjugate_bases[0][1], this->root_list[1]);                  // a*delta'^2 + b*delta'


  mul(this->conjugate_bases[1][0],to<PType>((this->defining_IBCF)[3]), this->root_list[2]);      // a*delta''
  this->conjugate_bases[1][1] = this->conjugate_bases[1][0];                              // a*delta''
  add(this->conjugate_bases[1][1],this->conjugate_bases[1][1], to<PType>( (this->defining_IBCF)[2]) ); // a*delta'' + b
  mul(this->conjugate_bases[1][1], this->conjugate_bases[1][1], this->root_list[2]);                    // a*delta''^2 + b*delta''
}

template <typename Type, typename PType>
void CubicOrderReal<Type, PType> :: get_real_value(PType & newVal, Type &U, Type &X, Type &Y, Type &D, int conj){

    if (this->is_real() && (conj != 0)){
       --conj;
       mul(newVal, to<PType>(X), this->get_conjugate_bases(conj,0));  // x*rho1'

       add(newVal, newVal, to<PType>(U));
       mul(this->order_temp, to<PType>(Y), this->get_conjugate_bases(conj,1)); // y* rho2'
       add(newVal, newVal, this->order_temp);
    }
    else{
      mul(newVal, to<PType>(X), this->get_rho1());
      add(newVal, newVal, to<PType>(U));
      mul(this->order_temp, to<PType>(Y), this->get_rho2());
      add(newVal, newVal, this->order_temp);
    }

    div(newVal, newVal, to<PType>(D));

}






/*
// This is a working version of the Voronoi algorithm for Fundamental units, but it
// has been moved to the BasicVornoi class

template <typename Type, typename PType>
void CubicOrderReal<Type, PType> :: compute_fundamental_unit(){

      // We use this->x_cycle and adj_min_vec to hold the cycle of lattices and corresponding minima

  //////////////////////////////////////////////////////////////////////////////
  //                          Local variables:
  //////////////////////////////////////////////////////////////////////////////
      //LatticeBasis L,Ltest;             // Variable for the iterated lattice bases
      //NMatrix epsilon1,                 // the fundamental unit to be.
      //          adj_minimum ,                // the minima adjacent to 1
      //          psi_bar;                //eventually holds second fund. unit.
      CubicElement<Type, PType> epsilon1 = CubicElement<Type, PType>(this);
      CubicElement<Type, PType> adj_minimum = CubicElement<Type, PType>(this);
      CubicElement<Type, PType> epsilon2 = CubicElement<Type, PType>(this);

      bool completeCycle;                // flag to indicate whether one full lattice
                                          // period has been traversed X and Z cycles resp.
      int cycleSize1 = 0;                 // counts the iterations in the XCycle

      // both of these are initialized as the identity, L1 is the matrix on which we
      // operate, L2 holds the 'adjacent' ideal, which gets swapped to L1 each loop
      CubicIdeal<Type, PType> L1 = CubicIdeal<Type, PType>(this);
      CubicIdeal<Type, PType> L2 = CubicIdeal<Type, PType>(this);
  ///////////////////////////////////////////////////////////////////////////////


      // initialize the variable which indicates when a full period has been traversed.
      completeCycle = false;

      ////////////////////////////////////////////////////////////////////////////
      //                                                                        //
      //                    Step 2: Voronoi along the X-axis                    //
      //                                                                        //
      ////////////////////////////////////////////////////////////////////////////
      while (!completeCycle){
          std::cout << "-----------------------------------------------------" << std::endl;
          std::cout << "----------------XCycle iteration: " << cycleSize1 << "------------------"<< std::endl;
          std::cout << "-----------------------------------------------------" << std::endl;
          ++cycleSize1;

          // make sure L1 is in canonical form, the push it onto x_cycle
          L1.make_canonical();
          this->x_cycle.push_back(L1);
          // convert to Voronoi basis
          L1.make_voronoi_basis();



          //divide out by adj_minimum and then swap so that L1 is the new lattice/ideal
          L1.adjacent_ideal(L2, adj_minimum);
          std::swap(L1,L2);
          // push adj_minima (which should be theta_g) onto adj_minima_vec
          adj_minima_vec.push_back(adj_minimum );

          // compare with the previous entries of x_cycle
          for (int i = 0; i < this->x_cycle.size(); i++){
              //std::cout << "Comparing with stored lattice :" << i << "..." << std::endl;
              completeCycle = ::is_equal(L1, this->x_cycle[i]);

              // Once a lattice L1 is found which matches a previous lattice L, we delete the lattices preceding L
              if (completeCycle){

                  std::cout << "Smallest match found with lattice " << i << std::endl;
                  std::cout << " Number of Lattices in this->x_cycle"<<this->x_cycle.size() << std::endl;

                  this->x_cycle.erase(this->x_cycle.begin(), this->x_cycle.begin()+i);
                  adj_minima_vec.erase(adj_minima_vec.begin(), adj_minima_vec.begin()+i);

                  std::cout << "Lattices remaining: " <<  this->x_cycle.size() << std::endl;

                  break;    // breaks out of the x_cycle comparison loop and enter step 4
              }
          } // end of x_cycle comparisons loop

      } //closes Step2 while loop
      std::cout << "-----------------------------------------------------" << std::endl;
      std::cout << "-------------------XCycle complete-------------------" << std::endl;
      std::cout << "-----------------------------------------------------" << std::endl;
      ////////////////////////////////////////////////////////////////////////////
      //                                                                        //
      //                    Step 4: Initialize epsilon2 and rotate axis          //
      //                                                                        //
      ////////////////////////////////////////////////////////////////////////////

      //At this point psi_bar (aka epsilon2) is equal to 1.
      // this should rotate the roots, and recompute the rho_i and conjugates
      this->roots_swap_position(0,2);


      // it should be that L1 is already equal to the 0th entry of x_cycle
      // the code below is unnecessary
      //for (int i = 0; i < 3; ++i){
      //  for (int j = 0; j < 2; ++j){
      //    L1.coefficientMatrix[i][j] = this->x_cycle[0].coefficientMatrix[i][j];
      //  }
      //} //end initializing loop for L
      //L1.mainDenominator = this->x_cycle[0].mainDenominator; //redundant, should always be 1

      #ifdef DEBUG
        std::cout << "\n" <<"FundamentalUnit: Lbar: " << std::endl;
        std::cout << L1.get_coeff(0,1) << "   " << L1.get_coeff(0,2) << " " << std::endl;
        std::cout << L1.get_coeff(1,1) << "   " << L1.get_coeff(1,2) << " " << std::endl;
        std::cout << L1.get_coeff(2,1) << "   " << L1.get_coeff(2,2) << " " << std::endl;
        std::cout << L1.get_coeff(0,0) << std::endl;
      #endif

      // reset the bool for the Z-cycle
      completeCycle = false;

      ////////////////////////////////////////////////////////////////////////////
      //                                                                        //
      //                    Step 5: Voronoi along the Z-axis                    //
      //                                                                        //
      ////////////////////////////////////////////////////////////////////////////
      while (!completeCycle){

          // convert to Voronoi Basis and extract theta_g into adj_min
          L1.make_voronoi_basis();

          // Obtain the adjacent ideal/lattice and swap so it is in L1.
          L1.divide_adjacent(L2, adj_minimum);
          std::swap(L1,L2);

          // update psi_bar (epsilon2)
          ::mul(epsilon2, epsilon2, adj_minimum);

          // compare with the lattices of x_cycle
          for (int i = 0; i < this->x_cycle.size(); i++){
              std::cout << "Comparing with stored lattice :" << i << "..." << std::endl;
              completeCycle = ::is_equal(L1, this->x_cycle[i]);


              // Once we've found a repeat lattice, we cave computed epsilon2
              if (completeCycle){
                  //std::cout << "epsilon2 " << epsilon2.get_u() << " "<<epsilon2.get_x() << " "<< epsilon2.get_y() << std::endl;

                  // Using the cycle size, we compute epsilon1 as the product of all the elements of the vector adj_minima_vec
                  // we also compute epsilon_2 using a psi and an intermediate product of the elements of adj_minima_vec,
                  // which is captured by the if statement below
                  for (int j = 0; j< this->x_cycle.size(); ++j){
                      if (j ==i){
                          //std::cout << "Intermediate theta: " << epsilon1.get_u() << " "<<epsilon1.get_x() << " "<< epsilon1.get_y() << std::endl;
                          epsilon1.inverse(adj_minimum);            //adj_minimum  = Inverse(epsilon1);
                          ::mul(epsilon2, epsilon2, adj_minimum);     //MatrixMultiplication(epsilon2, adj_minimum );
                      }
                      ::mul(epsilon1, epsilon1, adj_minima_vec[j]);
                  }
                  break;
              }
          }//end for loop

      }// end while loop


      //Undoes the original rotation.
      this->roots_swap_position(0,2);
      this->fundamentalUnits.push_back(epsilon1);
      this->fundamentalUnits.push_back(epsilon2);
      PType eps1, eps2;
      this->fundamentalUnits[0].get_real_value(eps1); this->fundamentalUnits[1].get_real_value(eps2);
      std::cout << "Epsilon_1 " << epsilon1.get_u() << " "<<epsilon1.get_x() << " "<< epsilon1.get_y() << " : " <<  eps1 << std::endl;
      std::cout << "Epsilon_2 " << epsilon2.get_u() << " "<<epsilon2.get_x() << " "<< epsilon2.get_y() << " : " <<  eps2 << std::endl;



}*/

































#endif
