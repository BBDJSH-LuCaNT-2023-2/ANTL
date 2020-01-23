#ifndef ANTL_VORONOI_REAL_CPP
#define ANTL_VORONOI_REAL_CPP

#include "../../include/ANTL/Cubic/VoronoiReal.hpp"
#include "../../include/ANTL/Cubic/CubicOrderReal.hpp"
template <typename Type, typename PType>

void VoronoiReal<Type, PType> :: make_voronoi_basis(CubicIdeal<Type, PType> & ideal1){


      this->real_ord = dynamic_cast<CubicOrderReal<Type, PType>*>(ideal1.get_order());
      //std::cout << "Main root "this->real_ord->get_root1() << std::endl;

      //std::cout << "rho1' and rho2' "<< this->real_ord->get_conjugate_bases(0,0) << " " << this->real_ord->get_conjugate_bases(0,1) << std::endl;
      //std::cout << "rho1' and rho2'' "<< this->real_ord->get_conjugate_bases(1,0) << " " << this->real_ord->get_conjugate_bases(1,1) << std::endl;
      // this->omegaMatrix[3][5];  //to hold the coordinates of omega_1 ... omega_5
                                   //which are the candidates for the 5 puncture thm

      int smallestGammaIndex;         // for use in step 11
      PType minGamma, GammaValue;     // for use in step 11


      // initialize omegamatrix whenever this function is called
      for (int i = 0; i < 5; ++i ){
        this->omegaDecision[i] = true;
      }
      //////////////////////////////////////////////////////////
      //                                                      //
      //     Step 1: make sure basis is in canonical form     //
      //                                                      //
      //////////////////////////////////////////////////////////

      //std::cout << "VoronoiBasis: Initial input matrix: " << std::endl;
      //std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      //std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  "   Denom: " << ideal1.coeff_matrix[0][0] << std::endl;
      //std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;

      ideal1.make_canonical();

      #ifdef DEBUG
      std::cout << "VoronoiBasis: Canonical form: " <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  "   Denom: " << ideal1.coeff_matrix[0][0] << std::endl;
      std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;


      //PType ppl[2][2];
      //ideal1.puncture_lattice(ppl);
      //std::cout << "canonical puncture" << std::endl;
      //std::cout << ppl[0][0] << " " << ppl[0][1] << std::endl;
      //std::cout << ppl[1][0] << " " << ppl[1][1] << std::endl;
      #endif

      //////////////////////////////////////////////////////////
      //                                                      //
      //        Step 2: get the puncture lattice of B         //
      //                                                      //
      //////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////
      //                                                             //
      // Step 3: convert B, PuncB to a prepared basis, {1, phi, psi} //
      //         As usual we leave out the first column              //
      //                                                             //
      /////////////////////////////////////////////////////////////////
      //PuncB = puncture(B);
      //becomePrepared(B, PuncB);
      ideal1.make_prepared();                 // accomplishes both Step 2 and 3

      #ifdef DEBUG
      std::cout << "VoronoiBasis: Prepared lattice: " <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;

      std::cout << "VoronoiBasis: puncture (Prepared): " << std::endl;
      std::cout <<  ideal1.p_lat[0][0] << "  " << ideal1.p_lat[0][1] <<  std::endl;
      std::cout <<  ideal1.p_lat[1][0] << "  " << ideal1.p_lat[1][1] <<  std::endl;
      #endif

      //////////////////////////////////////////////////////////
      //                                                      //
      //        Step 4: Modify phi and psi                    //
      //                                                      //
      //////////////////////////////////////////////////////////
      #ifdef DEBUG
      ideal1.get_order()->get_real_value(this->alpha2, ideal1.coeff_matrix[0][1],ideal1.coeff_matrix[1][1],ideal1.coeff_matrix[2][1], ideal1.coeff_matrix[0][0], 2);
      ideal1.get_order()->get_real_value(this->alpha1, ideal1.coeff_matrix[0][2],ideal1.coeff_matrix[1][2],ideal1.coeff_matrix[2][2], ideal1.coeff_matrix[0][0], 1);
      std::cout << ideal1.coeff_matrix[0][1] << " " << ideal1.coeff_matrix[1][1 ]<< " " << ideal1.coeff_matrix[2][1] << " " << std::endl;
      std::cout << "Phi'' : " << this->alpha2 << std::endl;
      std::cout << "Psi' : " << this->alpha1 << std::endl;
      #endif

      if(ideal1.p_lat[1][0] < 0){
          // phi = phi-floor(phi')-1
          // since these are all constant, we only modify the constant value of Phi
          // however in subtractig, we need to make sure to multiply by the denimonator:
          // [ -floor(phi') -1)] / s
          ideal1.get_order()->get_real_value(this->alpha2, ideal1.coeff_matrix[0][1],ideal1.coeff_matrix[1][1],ideal1.coeff_matrix[2][1], ideal1.coeff_matrix[0][0], 1);
          ideal1.coeff_matrix[0][1] = ideal1.coeff_matrix[0][1]
            - ideal1.coeff_matrix[0][0]*to<Type>(floor(this->alpha2)) - ideal1.coeff_matrix[0][0];

          ideal1.get_order()->get_real_value(this->alpha2, ideal1.coeff_matrix[0][2],ideal1.coeff_matrix[1][2],ideal1.coeff_matrix[2][2], ideal1.coeff_matrix[0][0], 2);
          ideal1.coeff_matrix[0][2] = ideal1.coeff_matrix[0][2]
            - ideal1.coeff_matrix[0][0]*to<Type>(floor(this->alpha2))- ideal1.coeff_matrix[0][0];
      }else{
        //std::cout << "VoronoiBasis: Step 4 option 2 " << std::endl;

        ideal1.get_order()->get_real_value(this->alpha2, ideal1.coeff_matrix[0][1],ideal1.coeff_matrix[1][1],ideal1.coeff_matrix[2][1], ideal1.coeff_matrix[0][0], 2);
        ideal1.coeff_matrix[0][1] = ideal1.coeff_matrix[0][1]
          - ideal1.coeff_matrix[0][0]*to<Type>(floor(this->alpha2)) - ideal1.coeff_matrix[0][0];

        ideal1.get_order()->get_real_value(this->alpha2, ideal1.coeff_matrix[0][2],ideal1.coeff_matrix[1][2],ideal1.coeff_matrix[2][2], ideal1.coeff_matrix[0][0], 1);
        ideal1.coeff_matrix[0][2] = ideal1.coeff_matrix[0][2]
          - ideal1.coeff_matrix[0][0]*to<Type>(floor(this->alpha2)) - ideal1.coeff_matrix[0][0];
      }
      //std::cout << "VoronoiBasis: After Step 4" << std::endl;
      //std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      //std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  std::endl;
      //std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;
      /////////////////////////////////////////////////////////////
      //                                                         //
      //        Step 5: Define the omega_i                       //
      //                                                         //
      // take care that they each carry around ideal1.coeff_matrix[0][0] //
      /////////////////////////////////////////////////////////////
      this->omegaMatrix[0][0] = ideal1.coeff_matrix[0][1];    //phi
      this->omegaMatrix[1][0] = ideal1.coeff_matrix[1][1];
      this->omegaMatrix[2][0] = ideal1.coeff_matrix[2][1];

      this->omegaMatrix[0][1] = ideal1.coeff_matrix[0][2];    //psi
      this->omegaMatrix[1][1] = ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][1] = ideal1.coeff_matrix[2][2];

      this->omegaMatrix[0][2] = ideal1.coeff_matrix[0][1]-ideal1.coeff_matrix[0][2];    //phi-psi
      this->omegaMatrix[1][2] = ideal1.coeff_matrix[1][1]-ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][2] = ideal1.coeff_matrix[2][1]-ideal1.coeff_matrix[2][2];

      this->omegaMatrix[0][3] = ideal1.coeff_matrix[0][1]+ideal1.coeff_matrix[0][2];    //phi+psi
      this->omegaMatrix[1][3] = ideal1.coeff_matrix[1][1]+ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][3] = ideal1.coeff_matrix[2][1]+ideal1.coeff_matrix[2][2];

      this->omegaMatrix[0][4] = (2*ideal1.coeff_matrix[0][1])+ideal1.coeff_matrix[0][2];  //2phi +psi
      this->omegaMatrix[1][4] = (2*ideal1.coeff_matrix[1][1])+ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][4] = (2*ideal1.coeff_matrix[2][1])+ideal1.coeff_matrix[2][2];



      if (ideal1.is_reduced() == char(1)){


      //////////////////////////////////////////////////////////
      //                                                      //
      //        Step 6: Narrow down the list of candidates    //
      //                                                      //
      //////////////////////////////////////////////////////////

      // Step 6 indicates an equality involving the quantity A of the Hessian.
      // Given by b^2-3ac
      sqr(this->alpha2, to<PType>(ideal1.get_order()->get_coeff(2))) ;
      mul(this->alpha1, to<PType>(ideal1.get_order()->get_coeff(3)), to<PType>(ideal1.get_order()->get_coeff(1) ));
      mul(this->alpha1, this->alpha1, to<PType>(3) );
      sub(this->alpha2, this->alpha2, this->alpha1);      //b^2-3ac

      sqr(this->alpha1, to<PType>(ideal1.coeff_matrix[0][0]) );
      mul(this->alpha1, this->alpha1, to<PType>(49) );
      div(this->alpha1, this->alpha1, to<PType>(4) );
      if((ideal1.p_lat[0][1] >to<PType>(1)) ||
        ( this->alpha2 > ((49/4.0)*this->alpha1)) || IsOne(ideal1.coeff_matrix[0][0]) ){
          this->omegaDecision[3] = false;
          this->omegaDecision[4] = false;
      }
      else if (ideal1.p_lat[0][1] > PType(0.5)){
          this->omegaDecision[4] = false;
      }
      #ifdef DEBUG
      std::cout << "VoronoiBasis: Step 6 remaining candidates: "<< this->omegaDecision[0] << this->omegaDecision[1] <<this->omegaDecision[2] << this->omegaDecision[3]<< this->omegaDecision[4] << std::endl;
      #endif
      //////////////////////////////////////////////////////////
      //                                                      //
      //        Step 7: Further eliminations                  //
      //                                                      //
      //////////////////////////////////////////////////////////
      abs(this->alpha1, ideal1.p_lat[1][1]);
      abs(this->alpha2, ideal1.p_lat[1][0]);
      if (this->alpha1 > (1 + 2*this->alpha2 ) ){
          this->omegaDecision[1] = false;
          this->omegaDecision[2] = false;
          this->omegaDecision[3] = false;
          this->omegaDecision[4] = false;
      }
      else if (this->alpha1 > (1 + this->alpha2 ) ){
          this->omegaDecision[1] = false;
          this->omegaDecision[2] = false;
          this->omegaDecision[3] = false;

      }
      else if (this->alpha1 > 1 ){
          this->omegaDecision[1] = false;
          this->omegaDecision[2] = false;
      }
      else if (this->alpha1 > (1 - this->alpha2 ) ){
          this->omegaDecision[2] = false;
          //std::cout << "T = 1,2,4,5" << "\n";
      }
      ///////////////////////////////////////////////////////////////
      //                                                           //
      //        Step 8: remove the invalidated candidates          //
      //        The method used already accomplished this step     //
      ///////////////////////////////////////////////////////////////

    } // this closes the stuff we skip if the ideal is not known to be reduced already.

      //////////////////////////////////////////////////////////////////
      //                                                              //
      //        Step 9: Removal of possibilities via                  //
      //                checking | floor(omega_i')-floor(omega_i'')|  //                                    //
      //////////////////////////////////////////////////////////////////
      for (int i = 1; i<5; ++i){

        ideal1.get_order()->get_real_value(this->alpha1,  this->omegaMatrix[0][i], this->omegaMatrix[1][i], this->omegaMatrix[2][i], ideal1.coeff_matrix[0][0], 1);
        ideal1.get_order()->get_real_value(this->alpha2, this->omegaMatrix[0][i], this->omegaMatrix[1][i], this->omegaMatrix[2][i], ideal1.coeff_matrix[0][0], 2);
        floor(this->alpha1, this->alpha1);
        floor(this->alpha2, this->alpha2);
        sub(this->alpha1, this->alpha1, this->alpha2);
        abs(this->alpha1, this->alpha1);
          if ( this->omegaDecision[i] && ( this->alpha1  > to<PType>(1) ))
             {
                this->omegaDecision[i] = false;
             }
      }
      #ifdef DEBUG
      std::cout << "VoronoiBasis: Step 9 remaining candidates: "<< this->omegaDecision[0] << this->omegaDecision[1] <<this->omegaDecision[2] << this->omegaDecision[3]<< this->omegaDecision[4] << std::endl;
      #endif
      //////////////////////////////////////////////////////////
      //                                                      //
      //        Step 10: Computing the gamma_i                //
      //                                                      //
      //////////////////////////////////////////////////////////

      if (ideal1.p_lat[1][0] < 0){
        p1 = 1; p2 = 2;
      }
      else if (ideal1.p_lat[1][0] > 0){
        p1 = 2; p2 = 1;
      }

          if (this->omegaDecision[2]){
                  ideal1.get_order()->get_real_value(this->alpha2,  this->omegaMatrix[0][2], this->omegaMatrix[1][2], this->omegaMatrix[2][2], ideal1.coeff_matrix[0][0], p1);
                  this->omegaMatrix[0][2] = this->omegaMatrix[0][2]
                      - ideal1.coeff_matrix[0][0] - to<Type>(floor(this->alpha2))*ideal1.coeff_matrix[0][0];
          }
          if (this->omegaDecision[3]){
                  ideal1.get_order()->get_real_value(this->alpha2,  this->omegaMatrix[0][3], this->omegaMatrix[1][3], this->omegaMatrix[2][3], ideal1.coeff_matrix[0][0], p2);
                  this->omegaMatrix[0][3] = this->omegaMatrix[0][3]
                      - ideal1.coeff_matrix[0][0] - to<Type>(floor(this->alpha2))*ideal1.coeff_matrix[0][0];
          }
          abs(this->alpha1,ideal1.p_lat[1][1] );
          if (this->omegaDecision[4] && (this->alpha1 > to<PType>(1))){
                  ideal1.get_order()->get_real_value(this->alpha2,  this->omegaMatrix[0][4], this->omegaMatrix[1][4], this->omegaMatrix[2][4], ideal1.coeff_matrix[0][0], p2);
                  this->omegaMatrix[0][4] = this->omegaMatrix[0][4]
                      - ideal1.coeff_matrix[0][0] - to<Type>(floor(this->alpha2))*ideal1.coeff_matrix[0][0];
          }

//      else if (ideal1.p_lat[1][0] > 0){
//        if (this->omegaDecision[2]){
//                this->omegaMatrix[0][2] = this->omegaMatrix[0][2]
//                    - ideal1.coeff_matrix[0][0] - to<Type>(floor(NumericalConjugate(2, this->omegaMatrix[0][2], this->omegaMatrix[1][2], this->omegaMatrix[2][2], ideal1.coeff_matrix[0][0])))*ideal1.coeff_matrix[0][0];
//        }
//        if (this->omegaDecision[3]){
//                this->omegaMatrix[0][3] = this->omegaMatrix[0][3]
//                    - ideal1.coeff_matrix[0][0] - to<Type>(floor(NumericalConjugate(1, this->omegaMatrix[0][3], this->omegaMatrix[1][3], this->omegaMatrix[2][3], ideal1.coeff_matrix[0][0])))*ideal1.coeff_matrix[0][0];
//        }
//        if (this->omegaDecision[4] && (ANTL::abs(ideal1.p_lat[1][1]) > 1)){
//                this->omegaMatrix[0][4] = this->omegaMatrix[0][4]
//                    - ideal1.coeff_matrix[0][0] - to<Type>(floor(NumericalConjugate(1, this->omegaMatrix[0][4], this->omegaMatrix[1][4], this->omegaMatrix[2][4], ideal1.coeff_matrix[0][0])))*ideal1.coeff_matrix[0][0];
//        }
//      }
      //std::cout <<  NumericalConjugate(0, this->omegaMatrix[0][0], this->omegaMatrix[1][0], this->omegaMatrix[2][0], ideal1.coeff_matrix[0][0]) << std::endl;
      //std::cout <<  NumericalConjugate(0, this->omegaMatrix[0][1], this->omegaMatrix[1][1], this->omegaMatrix[2][1], ideal1.coeff_matrix[0][0]) << std::endl;

      ///////////////////////////////////////////////////////////
      //                                                       //
      //        Step 11: Determine the smallest valid gamma_i  //
      //                                                       //
      ///////////////////////////////////////////////////////////
      smallestGammaIndex = 0;

      ideal1.get_order()->get_real_value(minGamma,  this->omegaMatrix[0][0], this->omegaMatrix[1][0], this->omegaMatrix[2][0], ideal1.coeff_matrix[0][0], 0);
      for (int i = 1; i < 5; ++i){
          if (this->omegaDecision[i]){
              ideal1.get_order()->get_real_value(GammaValue,  this->omegaMatrix[0][i], this->omegaMatrix[1][i], this->omegaMatrix[2][i], ideal1.coeff_matrix[0][0], 0);
              if (GammaValue < minGamma){
                  smallestGammaIndex = i;
                  minGamma = GammaValue;
              }
          }
      }
      if(ideal1.is_reduced() != char(1)){
        if(minGamma > to<PType>(1)){
          ideal1.set_reduction_state(1);
        } else{
          ideal1.set_reduction_state(0);
        }
      }
      #ifdef DEBUGVORONOI
      std::cout <<  "    VoronoiBasis: Index of the adjacent min:  " << smallestGammaIndex<< std::endl;
      #endif
      //////////////////////////////////////////////////////////
      //                                                      //
      //        Step 12: get the puncture lattice of B         //
      //                                                      //
      //////////////////////////////////////////////////////////
      ideal1.coeff_matrix[0][1] = this->omegaMatrix[0][smallestGammaIndex];
      ideal1.coeff_matrix[1][1] = this->omegaMatrix[1][smallestGammaIndex];
      ideal1.coeff_matrix[2][1] = this->omegaMatrix[2][smallestGammaIndex];

      if (smallestGammaIndex == 0){
          ideal1.coeff_matrix[0][2] = this->omegaMatrix[0][1];
          ideal1.coeff_matrix[1][2] = this->omegaMatrix[1][1];
          ideal1.coeff_matrix[2][2] = this->omegaMatrix[2][1];
      }else{
          ideal1.coeff_matrix[0][2] = this->omegaMatrix[0][0];
          ideal1.coeff_matrix[1][2] = this->omegaMatrix[1][0];
          ideal1.coeff_matrix[2][2] = this->omegaMatrix[2][0];
      }
      #ifdef DEBUGVORONOI
      std::cout << "\n"<< "    VoronoiBasis: Final Output: " << std::endl;
      std::cout <<  "    " << ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      std::cout <<  "    " <<ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  std::endl;
      std::cout <<  "    " <<ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;
      std::cout << "\n"<< std::endl;
      #endif

  }//close function





#endif // guard
