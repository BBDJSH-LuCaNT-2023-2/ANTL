#ifndef ANTL_BASIC_VORONOI_CPP
#define ANTL_BASIC_VORONOI_CPP

#include "../../../include/ANTL/Cubic/FundamentalUnits/BasicVoronoi.hpp"



template<typename Type, typename PType>
BasicVoronoi<Type, PType> :: BasicVoronoi(){

}



template<typename Type, typename PType>
void BasicVoronoi<Type, PType> :: fundamental_unit_complex(std::vector<CubicElement<Type, PType>> &unitvec, CubicOrder<Type, PType> * ord){


    // initialize matrix L0 as the identity
    CubicIdeal<Type, PType> L_identity = CubicIdeal<Type, PType>(ord);

    // ensure that L2 and L1 are the identity ideal of ord
    // ditto for epsilon and adj_minimum
    L1.set_order(ord); L2.set_order(ord);
    L1.make_identity(); L2.make_identity();
    epsilon1.set_order(ord); epsilon1.assign(Type(1));
    adj_minimum.set_order(ord); adj_minimum.assign(Type(1));

    completeCycle = false;

    int rounds=0;

    do {

      std::cout << " *************************************** " << std::endl;
      std::cout << " *************************************** " << std::endl;

          std::cout << "Fund.Unit Iteration Number --- " << rounds << std::endl;
          ++rounds; //prints and increases counter.


          // compute the VoronoiBasis (1, theta_g, theta_h),
          //  where theta_g is the minima adjacent to 1 in L
          //L.make_voronoi_basis();

  //std::cout << "FundamentalUnit2: After voronoiBasis: " << std::endl;
  //std::cout << L.coefficientMatrix[0][0] << "   " << L.coefficientMatrix[0][1] << " " << std::endl;
  //std::cout << L.coefficientMatrix[1][0] << "   " << L.coefficientMatrix[1][1] << " " << std::endl;
  //std::cout << L.coefficientMatrix[2][0] << "   " << L.coefficientMatrix[2][1] << " " << std::endl;
  //std::cout << L.mainDenominator << std::endl;
  //std::cout << " ****************************** " << std::endl;
  //std::cout << "FundamentalUnit2: Is VB the same? "<< std::endl;
  //        if (compareLattice3(Ltest,L)!= 1){
  //          std::cout << "ERROR: VoronoiBasis has returned a different lattice. ABORT";
  //          break;
  //        }
  //        else{
  //          std::cout << "VoronoiBasis has returned the same lattice, continuing:";
  //        }
  //std::cout << " ****************************** " << std::endl;


          // If proper, after this call, L2 is the adjacent lattice, and L1. spare_ideal_element
          // will contain the relative minima adjacent to 1.

          L1.adjacent_ideal(L2, adj_minimum );

          // mulitiply epsilon by the min adjacent to 1.
          ::mul(epsilon1, epsilon1, adj_minimum );

          std::swap(L1, L2);

          completeCycle = ::is_equal(L1, L_identity);

      } while( !completeCycle );

      std::cout << "Fundamental Unit computed: "<< std::endl;
      unitvec.push_back(epsilon1);


};

template<typename Type, typename PType>
void BasicVoronoi<Type, PType> :: fundamental_unit_real(std::vector<CubicElement<Type, PType>> &unitvec, CubicOrder<Type, PType> * ord){
        // We use this->x_cycle and adj_min_vec to hold the cycle of lattices and corresponding minima

    //////////////////////////////////////////////////////////////////////////////
    //                          Local variables:
    //////////////////////////////////////////////////////////////////////////////
          //set epsilon1, epsilon2, and adj min to the 1 element of ord
          epsilon1.set_order(ord); epsilon1.assign(Type(1));
          adj_minimum.set_order(ord); adj_minimum.assign(Type(1));
          epsilon2.set_order(ord); epsilon2.assign(Type(1));

          completeCycle = false;              // flag to indicate whether one full lattice
                                              // period has been traversed.

          int cycleSize1 = 0;                 // counts the iterations in the XCycle

          // ensure that L2 and L1 are the identity ideal of ord
          // ditto for epsilon and adj_minimum
          L1.set_order(ord); L2.set_order(ord);
          L1.make_identity(); L2.make_identity();
    ///////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //                    Step 2: Voronoi along the X-axis                    //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        while (!completeCycle){
            std::cout << "-----------------------------------------------------" << std::endl;
            std::cout << "----------------XCycle iteration: " << cycleSize1 << "------------------" << std::endl;
            //std::cout << "-----------------------------------------------------" << std::endl;
            ++cycleSize1;

            // make sure L1 is in canonical form, the push it onto x_cycle
            L1.make_canonical();
            this->x_cycle.push_back(L1);
            // convert to Voronoi basis
            //L1.make_voronoi_basis();



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

                    std::cout << " MATCH: Lattice " << i << std::endl;
                    std::cout << " Number of Lattices in this->x_cycle: "<<this->x_cycle.size() << std::endl;

                    this->x_cycle.erase(this->x_cycle.begin(), this->x_cycle.begin()+i);
                    adj_minima_vec.erase(adj_minima_vec.begin(), adj_minima_vec.begin()+i);

                    std::cout << " Erasing Preperiod: Lattices remaining: " <<  this->x_cycle.size() << std::endl;

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

        ord->roots_swap_position(0,2);

        #ifdef DEBUGVORONOI
        std::cout << "after swapping roots" << std::endl;
        std::cout << ord->get_root1() << " " << ord->get_root2() << " " << ord->get_root3() << std::endl;
        std::cout << ord->get_rho1()<<std::endl; //" " << ord->get_conjugate_bases(0,0) << " " << ord->get_conjugate_bases(1,0) << std::endl;
        std::cout << ord->get_rho2() << std::endl;//" " << ord->get_conjugate_bases(0,1) << " " << ord->get_conjugate_bases(1,1) << std::endl;
        #endif
        // it should be that L1 is already equal to the 0th entry of x_cycle


        #ifdef DEBUG
          std::cout << "\n" <<"FundamentalUnit: Lbar: " << std::endl;
          std::cout << L1.get_coeff(0,1) << "   " << L1.get_coeff(0,2) << " " << std::endl;
          std::cout << L1.get_coeff(1,1) << "   " << L1.get_coeff(1,2) << " " << std::endl;
          std::cout << L1.get_coeff(2,1) << "   " << L1.get_coeff(2,2) << " " << std::endl;
          std::cout << L1.get_coeff(0,0) << std::endl;
        #endif

        // reset the bool for the Z-cycle
        completeCycle = false;
        std::cout << "print 1" << std::endl;
        ////////////////////////////////////////////////////////////////////////////
        //                                                                        //
        //                    Step 5: Voronoi along the Z-axis                    //
        //                                                                        //
        ////////////////////////////////////////////////////////////////////////////
        while (!completeCycle){

            // convert to Voronoi Basis and extract theta_g into adj_min
            L1.make_voronoi_basis();
            std::cout << "print 2" << std::endl;
            std::cout << L2.toString() << std::endl;
            std::cout << adj_minimum.toString() << std::endl;

            // Obtain the adjacent ideal/lattice and swap so it is in L1.
            L1.divide_adjacent(L2, adj_minimum);
            std::cout << "print 3" << std::endl;
            std::swap(L1,L2);

            // update psi_bar (epsilon2)
            ::mul(epsilon2, epsilon2, adj_minimum);

            // compare with the lattices of x_cycle
            for (int i = 0; i < this->x_cycle.size(); i++){
                std::cout << "Comparing with stored lattice: " << i << " ..." << std::endl;
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
        ord->roots_swap_position(0,2);
        unitvec.push_back(epsilon1);
        unitvec.push_back(epsilon2);
        PType eps1, eps2;
        unitvec[0].get_real_value(eps1); unitvec[1].get_real_value(eps2);

        std::cout << "BV: eps1 = " << epsilon1.get_u() << " "<<epsilon1.get_x() << " "<< epsilon1.get_y() << ", approx: " <<  eps1 << "     eps2 = " << epsilon2.get_u() << " "<<epsilon2.get_x() << " "<< epsilon2.get_y() << ", approx: " <<  eps2 << std::endl;

};
#endif
