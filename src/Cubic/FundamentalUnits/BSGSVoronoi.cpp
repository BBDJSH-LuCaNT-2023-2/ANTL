#ifndef ANTL_BSGS_VORONOI_CPP
#define ANTL_BSGS_VORONOI_CPP

#include "../../../include/ANTL/Cubic/FundamentalUnits/BSGSVoronoi.hpp"





template<typename Type, typename PType>
BSGSVoronoi<Type, PType> :: BSGSVoronoi(){

}

template<typename Type, typename PType>
void BSGSVoronoi<Type, PType> :: initialize_giant_baby(CubicOrder<Type, PType> * ord){

  NTL::clear(distance);
  babysteps.clear();
  x_cycle.clear();
  //this vector will hold the products of adj min, from 0 .. i
  adj_minima_vec.clear();

  abs(realtemp, to<PType>(ord->get_discriminant()));
  pow(giantBound, realtemp, PType(0.25));
  ceil(giantBound, giantBound);

  std::cout << "giantBound: " << giantBound << " Discriminant: " << ord->get_discriminant() << std::endl;
  babysteps.reserve(conv<size_t>(giantBound)); // makes the size of babysteps about disc^1/4
  hash_function.set_modulus(ANTL::to<Type>(giantBound));
  bigFoot.set_order(ord);
  L1.set_order(ord); L2.set_order(ord);
  L1.make_identity(); L2.make_identity();


  epsilon1.set_order(ord); epsilon1.assign(Type(1));
  adj_minimum.set_order(ord); adj_minimum.assign(Type(1));
  collision = false;

}


template<typename Type, typename PType>
void BSGSVoronoi<Type, PType> :: giant_step(CubicIdeal<Type, PType> & currentIdeal, CubicElement<Type, PType> & currentMinimum, PType & logval ){
  std::cout << "BSGS: Current Ideal \n" << currentIdeal.toString() << std::endl;

  //multiply current ideal with giant_step ideal
  mul(currentIdeal, currentIdeal, bigFoot);

  currentMinimum.get_real_value(this->realtemp);

  #ifdef DEBUG
  std::cout << "the real value of the current min: "<<this->realtemp << std::endl;
  giant_min.get_real_value(this->realtemp);
  std::cout << "the real value of the giant min: "<<this->realtemp << std::endl;
  #endif

  //multiply current minimum with giant_step minimum
  mul(currentMinimum, currentMinimum, giant_min);

  // reduce the ideal and placei n L2
  currentIdeal.reduce(L2, adj_minimum); //adj_minimum is holding the element we divide by to get a reduced lattice

  adj_minimum.get_real_value(this->realtemp);
  std::cout << "giant step: adj_min: " << this->realtemp << std::endl;
  // update the current minima
  mul(currentMinimum, currentMinimum, adj_minimum);

  // swap the currentIdeal for its reduced one
  std::swap(currentIdeal, L2);


  #ifdef DEBUG
  currentMinimum.get_real_value(this->realtemp);
  std::cout << "giant step: new minima: " << this->realtemp << std::endl;
  #endif

  log(realtemp, this->realtemp);
  add(logval, logval, this->realtemp);
  add(logval, logval, giant_logarithm);

}

template<typename Type, typename PType>
void BSGSVoronoi<Type, PType> :: fundamental_unit_complex(std::vector<CubicElement<Type, PType>> &unitvec, CubicOrder<Type, PType> * ord){

    this->initialize_giant_baby(ord);
    std::cout << "giant baby initialized" << std::endl;
    // initialize matrix L0 as the identity
    CubicIdeal<Type, PType> L_identity = CubicIdeal<Type, PType>(ord);

    Type rounds(0);

    do {

      //store the norm of the current lattice
      L1.integral_norm(this->norm_holder);

      // add the norm and L1 to our hash table

      // babysteps will actually just store the index at which to find the
      // the lattice in x_cycle
      babysteps.emplace(std::make_pair(this->norm_holder, conv<std::size_t>(rounds) ) );

      x_cycle.push_back(L1);
      adj_minima_vec.push_back(epsilon1);
      ++rounds; //prints and increases counter.

      L1.adjacent_ideal(L2, adj_minimum );

      // mulitiply epsilon by the min adjacent to 1.
      ::mul(epsilon1, epsilon1, adj_minimum );

      if (L2.is_one()){
        unitvec.push_back(epsilon1);
        std::cout << "FundamentalUnits complete" << std::endl;
        epsilon1.get_real_value(this->distance);
        log(this->distance, this->distance);
        std::cout << "Reg. estimate " << this->distance << std::endl;
        return;
      }
      std::swap(L1, L2);

      //completeCycle = ::is_equal(L1, L_identity);
    } while( rounds <= to<Type>(giantBound) );



    // For testing purposes, prints the babystep lattices
    std::cout << "Baby steps computed: "<< std::endl;
    std::cout << "This is our baby step table: \n \n" << std::endl;
    for (int xc =0; xc < x_cycle.size(); xc++){
      std::cout << x_cycle[xc].toString() << std::endl;
    }


    // We assign the variables giant step and giant_min
    // then assign the logarithm value. We now search for a collision using giant steps


    this->bigFoot.assign(L1);
    std::cout << "Our bigFoot is the ideal: \n" << bigFoot.toString() << std::endl;

    this->giant_min.assign(epsilon1);


    giant_min.get_real_value(giant_logarithm);
    log(giant_logarithm, giant_logarithm);
    epsilon1.get_real_value(this->realtemp);
    log(realtemp, realtemp);
    add(this->distance, this->distance, realtemp);
      //unitvec.push_back(epsilon1);
      int gs = 1;
    collision = false;

    do {

      giant_step(L1, epsilon1, distance); //returns a lattice equivalent
      L1.integral_norm(norm_holder);

      if (babysteps.count(norm_holder) >= 1){

        auto possible_matches = babysteps.equal_range(norm_holder);

        for (auto it = possible_matches.first; it != possible_matches.second; ++it) {

          std::cout << it->second<< std::endl;
          std::cout << x_cycle[it->second].toString();
          if( is_equal(x_cycle[it->second], L1 )) {

            // More debug statements
            std::cout <<  "we found a collision!" << std::endl;
            std::cout << L1.toString() << std::endl;
            std::cout << x_cycle[it->second].toString() << std::endl;

            collision = true;

            // Since we may have overstepped, this division is like going
            // backwards to the identity
            div(epsilon1, epsilon1, adj_minima_vec[it->second]);
            adj_minima_vec[it->second].get_real_value(realtemp);
            log(realtemp,realtemp);
            sub(distance, distance, realtemp);

          }
        }
      }

      std::cout << gs << " giant steps taken thus far." << std::endl; gs++;

    } while( !collision);
unitvec.push_back(epsilon1);
std::cout << "FundamentalUnits complete" << std::endl;
std::cout << "Reg. estimate " << this->distance << std::endl;
};

template<typename Type, typename PType>
void BSGSVoronoi<Type, PType> :: fundamental_unit_real(std::vector<CubicElement<Type, PType>> &unitvec, CubicOrder<Type, PType> * ord){

std::cout<< "Warning, we haven't implemented this for totally real cubics! This is just basic Voronoi!" << std::endl;
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
        ord->roots_swap_position(0,2);

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
        ord->roots_swap_position(0,2);
        unitvec.push_back(epsilon1);
        unitvec.push_back(epsilon2);
        PType eps1, eps2;
        unitvec[0].get_real_value(eps1); unitvec[1].get_real_value(eps2);
        std::cout << "Epsilon_1 " << epsilon1.get_u() << " "<<epsilon1.get_x() << " "<< epsilon1.get_y() << " : " <<  eps1 << std::endl;
        std::cout << "Epsilon_2 " << epsilon2.get_u() << " "<<epsilon2.get_x() << " "<< epsilon2.get_y() << " : " <<  eps2 << std::endl;

};
#endif
