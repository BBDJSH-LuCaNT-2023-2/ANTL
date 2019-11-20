#ifndef ANTL_VORONOI_COMPLEX_CPP
#define ANTL_VORONOI_COMPLEX_CPP

#include "../../include/ANTL/Cubic/VoronoiComplex.hpp"

template <typename Type, typename PType>
void VoronoiComplex<Type, PType> :: make_voronoi_basis(CubicIdeal<Type, PType> & ideal1){

      //Type omegaMatrix[3][5];  //to hold the coordinates of omega_1 ... omega_5
                                      //which are the cadidates from the 5 puncture thm
      PType Fbar[3][5];               //Holds the transformed matrix from step 6
      PType gammaMatrix[3][3];        //This is Gamma_C' from step 6
      PType* candidate;
      int index, s;

      //these are used in step 5 as precomputed values
      Type Irho1, Irho2;


      //Needed so that a column of omegaMatrix (which is an element of K)
      //can be interpreted as such. Needed for realEmbedding and cubicTrace
      // ZZVector3 elementHolder;  replaced by this->placeholder

      //this bool array will be to decide the correct omega
      for (int i = 0; i < 5; ++i ){
        this->omegaDecision[i] = true;
      }

      // initialization of the matrix Gamma, see pg 392 of CFG
      NTL::set(gammaMatrix[0][0]);
      NTL::clear(gammaMatrix[1][0]);
      NTL::set(gammaMatrix[2][0]);

      // Note that entries [1][1] and [1][2] should have 'i' in their expression,.
      // This is circumvented by observing delta'-delta'' = 2*Im(delta')*i
      gammaMatrix[0][1] =  ideal1.get_order()->get_rho1();                                             // rho1

      mul(gammaMatrix[1][1], to<PType>(ideal1.get_order()->get_coeff(3)), ideal1.get_order()->get_root3() );  // a*Im(delta')
      //gammaMatrix[1][1] = ( to<PType>(ideal1.get_order()->get_coeff(3))*ideal1.get_order()->get_root3() );    // a*Im(delta')

      NTL::clear(gammaMatrix[2][1]);
      sub(gammaMatrix[2][1],gammaMatrix[2][1], ideal1.get_order()->get_rho1());
      sub(gammaMatrix[2][1],gammaMatrix[2][1], to<PType>(ideal1.get_order()->get_coeff(2)));

      div(gammaMatrix[2][1], gammaMatrix[2][1], to<PType>(2));
      //gammaMatrix[2][1] = -(ideal1.get_order()->get_rho1() + to<PType>(ideal1.get_order()->get_coeff(2)))/2; // -(rho1 + b)/2
      gammaMatrix[0][2] =  ideal1.get_order()->get_rho2();                                             // rho2

      NTL::clear(gammaMatrix[1][2]);                                                                  // 0
      sub(gammaMatrix[1][2], gammaMatrix[1][2], to<PType>(ideal1.get_order()->get_coeff(3)));             // -a
      mul(gammaMatrix[1][2], gammaMatrix[1][2], ideal1.get_order()->get_root1());                     // -a*delta
      mul(gammaMatrix[1][2], gammaMatrix[1][2], ideal1.get_order()->get_root3());                     // -a*delta* im(delta')
      //gammaMatrix[1][2] = -(to<PType>(ideal1.get_order()->get_coeff(3))*(ideal1.get_order()->get_root1())*(ideal1.get_order()->get_root3()) );

      NTL::clear(gammaMatrix[2][2]);
      sub(gammaMatrix[2][2], gammaMatrix[2][2], to<PType>(ideal1.get_order()->get_coeff(1)));             // -c
      mul(gammaMatrix[2][2],gammaMatrix[2][2],2);                                                     //-2c
      sub(gammaMatrix[2][2],gammaMatrix[2][2], ideal1.get_order()->get_rho2());                       // -rho2 -2c
      div(gammaMatrix[2][2],gammaMatrix[2][2],2);                                                     // -(rho2 + 2c)/2
      //gammaMatrix[2][2] = -(ideal1.get_order()->get_rho2() + 2*to<PType>(ideal1.get_order()->get_coeff(1)))/2;  // -(rho2 + 2c)/2

      //std::cout << "gammaMatrix: " << std::endl;
      //std::cout <<  gammaMatrix[0][0] << "  " << gammaMatrix[0][1] << "  " << gammaMatrix[0][2]<< std::endl;
      //std::cout <<  gammaMatrix[1][0] << "  " << gammaMatrix[1][1] << "  " << gammaMatrix[1][2]<< std::endl;
      //std::cout <<  gammaMatrix[2][0] << "  " << gammaMatrix[2][1] << "  " <<gammaMatrix[2][2]<< std::endl;

      // Step 1: make sure basis is in canonical form
      std::cout << "VoronoiBasis: Initial input matrix: " << std::endl;
      std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;

      ideal1.make_canonical();
      std::cout << "VoronoiBasis: The canonical form: " <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;
      //std::cout << "VoronoiBasis: Check if canonical basis is the same lattice  "  << compareLattice3(B, Ltest)<<std::endl;


      // Step 2: Obtain the puncture lattice of B

      // Note that the CubicIdeal.make_prepared() function computes and stores
      // the corresponding puncture lattice into CubicIdeal.p_lat
      ideal1.make_prepared();

      // ideal1.puncture_lattice(temp_pb);
      // becomePrepared(ideal1, temp_pb);

      // Step 3: convert B, PuncB to a prepared basis, {1, phi, psi}
      // this step is performed above

      std::cout << "VoronoiBasis: Now Prepared: " <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  std::endl;
      std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;

      //std::cout << "VoronoiBasis: Check if prepared basis is the same lattice: "  << compareLattice3(B, Ltest)<<std::endl;
      //std::cout << "VoronoiBasis: Redundancy check (ensure sameness after comparison)" << std::endl;
      //std::cout <<  ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] <<  std::endl;
      //std::cout <<  ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] <<  std::endl;
      //std::cout <<  ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] <<  std::endl;
      //std::cout <<  ideal1.coeff_matrix[0][0] << std::endl;

      std::cout << "VoronoiBasis: PB's puncture Lattice:" << std::endl;
      std::cout <<  ideal1.p_lat[0][0] << "  " << ideal1.p_lat[0][1] <<  std::endl;
      std::cout <<  ideal1.p_lat[1][0] << "  " << ideal1.p_lat[1][1] <<  std::endl;
      #ifdef DEBUG
      std::cout << ( to<PType>(ideal1.coeff_matrix[0][1])
                    + to<PType>(ideal1.coeff_matrix[1][1]) *ideal1.get_order()->get_rho1()
                    + to<PType>(ideal1.coeff_matrix[2][1])*ideal1.get_order()->get_rho2())/to<PType>(ideal1.coeff_matrix[0][0]) << std::endl;
      std::cout << ( to<PType>(ideal1.coeff_matrix[0][2])
                    + to<PType>(ideal1.coeff_matrix[1][2]) *(ideal1.get_order()->get_rho1())
                    + to<PType>(ideal1.coeff_matrix[2][2])*ideal1.get_order()->get_rho2())/to<PType>(ideal1.coeff_matrix[0][0]) << std::endl;
      #endif

      // Step 4, formation ofto<Type>( the matrix F whose columns are basis reps of
      //  phi, psi, phi-psi, phi+psi, 2phi+psi

      //  These are referred to as omega_i,
      //  take care that they each carry around ideal1.coeff_matrix[0][0]


      this->omegaMatrix[0][0] = ideal1.coeff_matrix[0][1];    //phi
      this->omegaMatrix[1][0] = ideal1.coeff_matrix[1][1];
      this->omegaMatrix[2][0] = ideal1.coeff_matrix[2][1];

      this->omegaMatrix[0][1] = ideal1.coeff_matrix[0][2];    //psi
      this->omegaMatrix[1][1] = ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][1] = ideal1.coeff_matrix[2][2];

      this->omegaMatrix[0][2] = ideal1.coeff_matrix[0][1] - ideal1.coeff_matrix[0][2];    //phi-psi
      this->omegaMatrix[1][2] = ideal1.coeff_matrix[1][1] - ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][2] = ideal1.coeff_matrix[2][1] - ideal1.coeff_matrix[2][2];

      this->omegaMatrix[0][3] = ideal1.coeff_matrix[0][1]+ideal1.coeff_matrix[0][2];    //phi+psi
      this->omegaMatrix[1][3] = ideal1.coeff_matrix[1][1]+ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][3] = ideal1.coeff_matrix[2][1]+ideal1.coeff_matrix[2][2];

      this->omegaMatrix[0][4] = (2*ideal1.coeff_matrix[0][1])+ideal1.coeff_matrix[0][2];  //2phi +psi
      this->omegaMatrix[1][4] = (2*ideal1.coeff_matrix[1][1])+ideal1.coeff_matrix[1][2];
      this->omegaMatrix[2][4] = (2*ideal1.coeff_matrix[2][1])+ideal1.coeff_matrix[2][2];

      #ifdef DEBUG
      std::cout << "VoronoiBasis: Omega matrix" << std::endl;
      std::cout << this->omegaMatrix[0][0] << " " << this->omegaMatrix[0][1] << " " << this->omegaMatrix[0][2] << " " << this->omegaMatrix[0][3] << " " << this->omegaMatrix[0][4] << std::endl;
      std::cout << this->omegaMatrix[1][0] << " " << this->omegaMatrix[1][1] << " " << this->omegaMatrix[1][2] << " " << this->omegaMatrix[1][3] << " " << this->omegaMatrix[1][4] << std::endl;
      std::cout << this->omegaMatrix[2][0] << " " << this->omegaMatrix[2][1] << " " << this->omegaMatrix[2][2] << " " << this->omegaMatrix[2][3] << " " << this->omegaMatrix[2][4] << std::endl;
      #endif
      ///////////////////////////////////////////////////////////////////////////
      //
      //                                Step 5:
      //
      // Here we want to obtain the 'bar' versions of the omega_i,
      // This require the computation of [-zeta], and can be done only
      // with rational arithmetic.
      //
      ///////////////////////////////////////////////////////////////////////////

      // This is the method suggested by Hambleton and Williams, found on page 388
      // I don't actually believe I've implemented the condition correctly.
      abs(this->alpha2, to<PType>(ideal1.get_order()->get_discriminant()) );
      div(this->alpha1, to<PType>(4), to<PType>(3));
      pow(this->alpha1, this->alpha2, this->alpha1);

          if (this->alpha2 > to<PType>(16384) ){

              this->pm = to<Type>(ceil(to<PType>(213)*this->alpha1 ));
              //this->pm = to<Type>(ceil(to<PType>(213)*power(to<PType>(ideal1.get_order()->get_discriminant()),4/3) ));  // pm is playing the role of I in the text
              Irho1 = to<Type>(floor(to<PType>(this->pm) * ideal1.get_order()->get_rho1()));
              Irho2 = to<Type>(floor(to<PType>(this->pm)* ideal1.get_order()->get_rho2()));

              for (int i = 0; i <5; ++i){
                // if omega_i = (u + x*rho1 + y*rho2)/sigma
                // This formula is (I*(sigma -q) + x*Irho1 +y*Irho2) / (sigma*I*2)
                // where q = 2u - bx -2cy
                this->omegaMatrix[0][i] +=
                  // should floor automatically if all types are ZZ
                   ideal1.coeff_matrix[0][0]*(
                  ( this->pm*(ideal1.coeff_matrix[0][0]
                    - (2*this->omegaMatrix[0][i] - ideal1.get_order()->get_coeff(2)*this->omegaMatrix[1][i]
                    - 2*ideal1.get_order()->get_coeff(1)*this->omegaMatrix[2][i]) )          // I*(sigma - q)
                  + this->omegaMatrix[1][i]*Irho1
                  + this->omegaMatrix[2][i]*Irho2 )
                  /(2*ideal1.coeff_matrix[0][0]*this->pm));
              }
          }
          else{

            for (int i = 0; i < 5; ++i){
              // There are two other formulas for obtaining [-zeta ]
              // We use the formula 0.5*(Trace(theta) - theta) for zeta

              // Note: The elements of omegaMatrix
              // should all be divided by ideal1.coeff_matrix[0][0], we calculate what the
              // true value of [-zeta] is, but in order to properly add it to the ,
              // matrix elements, we have to put it over a common denominator.
              this->placeholder.set_order(ideal1.get_order());
              this->placeholder.assign(this->omegaMatrix[0][i],this->omegaMatrix[1][i],this->omegaMatrix[2][i], ideal1.coeff_matrix[0][0]);



              //std::cout << (0.5*( 3*realEmbedding(elementHolder,ideal1.coeff_matrix[0][0])
              //              -cubicTrace(elementHolder, ideal1.coeff_matrix[0][0]))
              //              -realEmbedding(elementHolder,ideal1.coeff_matrix[0][0]) )*to<PType>(ideal1.coeff_matrix[0][0])
              //              << std::endl;
              //QQ<Type> myrational;

              this->placeholder.trace(ideal1.rational_temp);                            // store the trace in dummy1

              div(this->alpha2, to<PType>(ideal1.rational_temp.getN()), to<PType>(ideal1.rational_temp.getD()) );

              this->placeholder.get_real_value(this->alpha0);


              //std::cout <<  "Value of -zeta" << i << "     " << (-0.5 *(this->alpha2
              //  -this->alpha0)) << std::endl;
              //std::cout << "Nearest integer (multiplied by denominator):      " << round (-0.5 *(this->alpha2
              //    -this->alpha0))*to<PType>(ideal1.coeff_matrix[0][0]) << std::endl;

              this->omegaMatrix[0][i] += to<Type>(
                round ( -0.5 *(this->alpha2
                      -this->alpha0))*to<PType>(ideal1.coeff_matrix[0][0]) );


            }
          }
          #ifdef DEBUG
          std::cout << "omegabar matrix" << std::endl;
          std::cout << this->omegaMatrix[0][0] << "  " << this->omegaMatrix[0][1] << "  " << this->omegaMatrix[0][2] << "  " << this->omegaMatrix[0][3] << "  " << this->omegaMatrix[0][4] << std::endl;
          std::cout << this->omegaMatrix[1][0] << "  " << this->omegaMatrix[1][1] << "  " << this->omegaMatrix[1][2] << "  " << this->omegaMatrix[1][3] << "  " << this->omegaMatrix[1][4] << std::endl;
          std::cout << this->omegaMatrix[2][0] << "  " << this->omegaMatrix[2][1] << "  " << this->omegaMatrix[2][2] << "  " << this->omegaMatrix[2][3] << "  " << this->omegaMatrix[2][4] << std::endl;
          #endif


  //        std::cout << "Output the arithemtic matrices of the candidates:" << std::endl;
  //        NMatrix candidateArith;
  //        LargeNumber fakedenom = to<Type>(1);
  //        LargeNumber secondOmega;
  //        for (int i = 0; i < 5; ++i){
  //            candidateArith = arithmeticMatrix(omegaMatrix[0][i],omegaMatrix[1][i],omegaMatrix[2][i], fakedenom);
  //            std::cout << candidateArith.mat[0][0]<< " " << candidateArith.mat[0][1] << " " <<  candidateArith.mat[0][2]<< " " << std::endl;
  //            std::cout << candidateArith.mat[1][0]<< " " << candidateArith.mat[1][1] << " " <<  candidateArith.mat[1][2]<< " " << std::endl;
  //            std::cout << candidateArith.mat[2][0]<< " " << candidateArith.mat[2][1] << " " <<  candidateArith.mat[2][2]<< " " << std::endl;
  //            std::cout << "-------------------" << std::endl;
  //            secondOmega = omegaMatrix[0][i]-to<Type>(1);
  //            candidateArith = arithmeticMatrix(secondOmega,omegaMatrix[1][i],omegaMatrix[2][i], fakedenom);
  //            std::cout << candidateArith.mat[0][0]<< " " << candidateArith.mat[0][1] << " " <<  candidateArith.mat[0][2]<< " " << std::endl;
  //            std::cout << candidateArith.mat[1][0]<< " " << candidateArith.mat[1][1] << " " <<  candidateArith.mat[1][2]<< " " << std::endl;
  //            std::cout << candidateArith.mat[2][0]<< " " << candidateArith.mat[2][1] << " " <<  candidateArith.mat[2][2]<< " " << std::endl;
  //            std::cout << "-------------------" << std::endl;
  //            std::cout << "-------------------" << std::endl;
  //        }


          // alternatively, we can use the formula 1/3 (xi_omega - Trace(omega)
          // The issue here is that xi can be obtained from our puncture lattice,
          // but it is specific for each i

//          omegaMatrix[0][0] += round(
//            ( puncLat[0][0] -
//            cubicTrace({omegaMatrix[0][0],omegaMatrix[1][0],omegaMatrix[2][0] },ideal1.coeff_matrix[0][0]))/3
//            );
//          omegaMatrix[0][1] += round(
//            ( puncLat[0][1] -
//            cubicTrace({omegaMatrix[0][1],omegaMatrix[1][1],omegaMatrix[2][1] },ideal1.coeff_matrix[0][0]))/3
//            );
//          omegaMatrix[0][2] += round(
//            ( puncLat[0][0] - puncLat[0][1] -
//            cubicTrace({omegaMatrix[0][2],omegaMatrix[1][2],omegaMatrix[2][2] },ideal1.coeff_matrix[0][0]))/3
//            );
//          omegaMatrix[0][3] += round(
//            ( puncLat[0][0] + puncLat[0][1]-
//            cubicTrace({omegaMatrix[0][3],omegaMatrix[1][3],omegaMatrix[2][3] },ideal1.coeff_matrix[0][0]))/3
//            );
//          omegaMatrix[0][4] += round(
//            ( 2*puncLat[0][0] + puncLat[0][1] -
//            cubicTrace({omegaMatrix[0][4],omegaMatrix[1][4],omegaMatrix[2][4] },ideal1.coeff_matrix[0][0]))/3
//            );


        //////////////////////////////////////////////////////////////////////////
        //
        //                                Step 6:
        //
        // Create the matrix Fbar
        // Note should also consider each of these entries as being divided
        // by B. mainDenominator
        //
        //////////////////////////////////////////////////////////////////////////
  //std::cout << "Fbar Matrix: " << std::endl;
          for (int j = 0 ; j <5; ++j){
              for(int i = 0; i < 3; ++i){
                  Fbar[i][j] = (gammaMatrix[i][0]*to<PType>(this->omegaMatrix[0][j])
                              + gammaMatrix[i][1]*to<PType>(this->omegaMatrix[1][j])
                              + gammaMatrix[i][2]*to<PType>(this->omegaMatrix[2][j]))/to<PType>(ideal1.coeff_matrix[0][0]);

              }
          }

  #ifdef DEBUG
  std::cout << Fbar[0][0] << " " << Fbar[0][1] << " " << Fbar[0][2] << " " << Fbar[0][3]<< " " << Fbar[0][4] << std::endl;
  std::cout << Fbar[1][0] << " " << Fbar[1][1] << " " << Fbar[1][2] << " " << Fbar[1][3]<< " " << Fbar[1][4] << std::endl;
  std::cout << Fbar[2][0] << " " << Fbar[2][1] << " " << Fbar[2][2] << " " << Fbar[2][3]<< " " << Fbar[2][4] << std::endl;
  std::cout << "log values " << std::endl;
  for (int k =0; k < 5; k++){
    if (Fbar[0][k] > 0)
      log(this->alpha1, Fbar[0][k]);
      std::cout << k << " " << this->alpha1 << std::endl;
  }


          std::cout << "t = zeta^2 + eta^2" << std::endl;
          for (int k = 0; k < 5; ++k){
            //std::cout << "entries: " << Fbar[1][k] << "  "<< Fbar[2][k]-1 << std::endl;
            //std::cout << "squares: " << power(Fbar[1][k],2) << "  "<< power(Fbar[2][k]-1,2) << std::endl;

            sub(this->alpha1,Fbar[2][k], 1 );               // (Fbar[2][k]-1))
            NTL::sqr(this->alpha1,this->alpha1);            // (Fbar[2][k]-1))^2
            NTL::sqr(this->alpha2, Fbar[1][k]);             // (Fbar[1][k])^2
            add(this->alpha1, this->alpha1, this->alpha2);  // (Fbar[1][k])^2 + (Fbar[2][k]-1))^2

            //PType testvalue = power(Fbar[1][k],2) + power((Fbar[2][k]-1),2);

            std::cout << k << " minus: " << this->alpha1 << ",  original: " << (this->alpha1  + 2*Fbar[2][k] -1) << std::endl;

          }
          #endif

        //////////////////////////////////////////////////////////////////////////
        //                                                                      //
        //                                IMPORTANT!!                           //
        //                                                                      //
        // Procedure branches depending on whether the initial ideal is reduced //
        // This should be indicated by the parameter 'reduced'                  //
        //                                                                      //
        //////////////////////////////////////////////////////////////////////////
        if (ideal1.is_reduced() == char(1)){

        //////////////////////////////////////////////////////////////////////////
        //
        //                                Step 7:
        //
        // Elimination of puncture candidates
        //
        //////////////////////////////////////////////////////////////////////////
          div(this->alpha2, to<PType>(3), to<PType>(2));
          SqrRoot(this->alpha2, this->alpha2);
          if (ideal1.p_lat[1][1] < this->alpha2 ){
              this->omegaDecision[4]= false;
              this->omegaDecision[3]= false;
              std::cout << "We've ruled out the last two candidates" << std::endl;
          }
          else if (ideal1.p_lat[0][1] > (1 - ideal1.p_lat[0][0])){
              this->omegaDecision[4]= false;
              std::cout << "We've ruled out the last candidate" << std::endl;

          }


        //////////////////////////////////////////////////////////////////////////
        //
        //                                Step 8,9:
        //
        // Step 8, 9 also eliminate 5 puncture candidates. Requires computation
        // of minimum of a cubic form, which is a heavy computation.
        // Partially implements, We are considering the minimum of forms.
        //
        //////////////////////////////////////////////////////////////////////////
        if (IsOne(ideal1.coeff_matrix[0][0])){
          this->omegaDecision[4]= false;
          this->omegaDecision[3]= false;
        }


      } // If the ideal is not reduced, I believe we cannot use these checks,
        // Everything else is the same since we just find the minimal element whose
        // Norm sits inside of N(1).
        // i.e we never check if the abs value of the element is smaller than 1 or not


        //////////////////////////////////////////////////////////////////////////
        //
        //                                Step 10:
        //
        // Check the radius of the remaining candidates w-1 and w
        // in order of minimal first coordinate. Find the smallest one
        //
        //
        //////////////////////////////////////////////////////////////////////////
          do{
              candidate = NULL;
              index = -1;

              //This for loop finds the minimal f_1i out of the remaining candidates
              for (int i = 0; i <5; ++i){
                  //if omega_i is still in the running, and the loop is just starting,
                  // select the first valid omega_i,
                  if ((index == -1) && (this->omegaDecision[i] == true)){
                      candidate = &Fbar[0][i];
                      index = i;
                  }
                  // If we've selected an omega, this loop goes through the remaining
                  // and yields the one with minimal \bar(f_{1j})
                  if (this->omegaDecision[i] == true && Fbar[0][i] < *candidate){
                      candidate = &Fbar[0][i];
                      index = i;
                  }

                  //std::cout << (*candidate) << std::endl;

              }

              // At this point, index denotes the omega with minimal real value.
              std::cout << "Smallest candidate is " << index << std::endl;
              this->omegaDecision[index] = false;     //10.b 'removal' of the index

              // Reusing the 0,0 entry of gammaMatrix
              //to store p(omega_i-1)
              // This tells us whether (omega_index) -1 lies in the cylinder C

              sub(gammaMatrix[0][0], Fbar[2][index], 1);                 // (Fbar[2][index]-1))
              NTL::sqr(gammaMatrix[0][0],gammaMatrix[0][0]);             // (Fbar[2][index]-1))^2
              NTL::sqr(this->alpha2, Fbar[1][index] );                   // (Fbar[1][index])^2
              add(gammaMatrix[0][0], gammaMatrix[0][0], this->alpha2);   // (Fbar[1][index])^2 + (Fbar[2][index]-1))^2
              //gammaMatrix[0][0] = power(Fbar[1][index],2) + power((Fbar[2][index]-1),2);

              // Part 10.c.


              // If (omega_index)-1 lies in C, this is the adjacent min
              if ( gammaMatrix[0][0] < to<PType>(1) ){
                  s = 1;   //this will be value 's' from step 10.c

                  break; //exit the do while, we've found the adjacent min.
              }
              // if the above is false, but omega_index lies in C, then there's some work to do.
              else if ( (gammaMatrix[0][0] + 2*Fbar[2][index] -1) < to<PType>(1)) {
                  s = 0;

                  // For the remaining indices j, check if the remaining omega_j-1
                  //are smaller than omega_index

                  for(int j = 0; j < 5; j++) {

                      if ((this->omegaDecision[j]) && ((Fbar[0][j]-1) < (Fbar[0][index]-s)) ){        // determine if the candidate is still in the running

                          sub(gammaMatrix[0][0], Fbar[2][j], 1);                     // (Fbar[2][index]-1))
                          NTL::sqr(gammaMatrix[0][0],gammaMatrix[0][0]);             // (Fbar[2][index]-1))^2
                          NTL::sqr(this->alpha2, Fbar[1][j] );                       // (Fbar[1][index])^2
                          add(gammaMatrix[0][0], gammaMatrix[0][0], this->alpha2);   // (Fbar[1][index])^2 + (Fbar[2][index]-1))^2

                          //gammaMatrix[0][0] = power(Fbar[1][j],2) + power((Fbar[2][j]-1),2);  // accomplished by the above 4 ops
                              // if smaller, check if in C.
                          if (gammaMatrix[0][0] < 1){
                              index = j; // j becomes new index
                              s = 1;     // tracks we are now dealing with omega_index-1
                          }
                      }
                  } // end for loop

                  break; // break out the while loop, we have found the adj. min
              }//close else if

          }
          while(index != -1);


          // Additional step when not reduced, check if the candidate is greater than 1
          // If not, then we are not reduced, but if it is, then we are reduced
          if(ideal1.is_reduced() != char(1)){
            sub(this->alpha2, Fbar[0][index], to<PType>(s));
            if(this->alpha2 > to<PType>(1)){
              ideal1.set_reduction_state(1);
            }
            else{
              ideal1.set_reduction_state(0);
            }
          }
        //////////////////////////////////////////////////////////////////////////
        //
        //                                Step 11:
        //
        // Once the adjacent minima is found, it is straightforward to pick the
        // final basis element. It's either bar(omega_1), or if we picked that already
        // bar(omega_2)
        //
        //////////////////////////////////////////////////////////////////////////


          if (index == 0){
              // Note: It is important to multiply s by ideal1.coeff_matrix[0][0],
              // We are 'adding 1' to the element, and without it, we might be adding (s/ideal1.coeff_matrix[0][0]) by accident
              ideal1.coeff_matrix[0][1] = this->omegaMatrix[0][0]-(s*ideal1.coeff_matrix[0][0]);
              ideal1.coeff_matrix[1][1] = this->omegaMatrix[1][0];
              ideal1.coeff_matrix[2][1] = this->omegaMatrix[2][0];
              ideal1.coeff_matrix[0][2] = this->omegaMatrix[0][1];
              ideal1.coeff_matrix[1][2] = this->omegaMatrix[1][1];
              ideal1.coeff_matrix[2][2] = this->omegaMatrix[2][1];
          }
          else{
              ideal1.coeff_matrix[0][1] = this->omegaMatrix[0][index]-(s*ideal1.coeff_matrix[0][0]);
              ideal1.coeff_matrix[1][1] = this->omegaMatrix[1][index];
              ideal1.coeff_matrix[2][1] = this->omegaMatrix[2][index];
              ideal1.coeff_matrix[0][2] = this->omegaMatrix[0][0];
              ideal1.coeff_matrix[1][2] = this->omegaMatrix[1][0];
              ideal1.coeff_matrix[2][2] = this->omegaMatrix[2][0];
          }

          std::cout << "VB New matrix (not canonical)" << std::endl;
          std::cout << ideal1.coeff_matrix[0][0] << "  "<< ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] << std::endl;
          std::cout << ideal1.coeff_matrix[1][0] << "  "<< ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] << std::endl;
          std::cout << ideal1.coeff_matrix[2][0] << "  "<< ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] << std::endl;

          ideal1.normalize();

std::cout << "***************** VoronoiEnd ******************"<< std::endl;
std::cout << "***************** ********** ******************"<< std::endl;
}



#endif // guard
