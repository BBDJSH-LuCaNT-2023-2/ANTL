#ifndef ANTL_VORONOI_METHODS_CPP
#define ANTL_VORONOI_METHODS_CPP

#include "../../include/ANTL/Cubic/VoronoiMethods.hpp"


template <typename Type, typename PType>
void VoronoiMethods<Type, PType>::make_prepared(CubicIdeal<Type, PType> & ideal1, PType plattice[2][2]){

      // Step 1:
      // Check if xi_mu is negative, flip the first row if it is
      if (plattice[0][0] < 0){

          plattice[0][0] = -plattice[0][0];
          plattice[1][0] = -plattice[1][0];

          ideal1.flip_basis_element(1);

      }

      // Check if xi_nu is negative, flip the second row if it is
      if (plattice[0][1] < 0){

          plattice[0][1] = -plattice[0][1];
          plattice[1][1] = -plattice[1][1];

          ideal1.flip_basis_element(2);

      }

      // STEP 2
      if ( (ANTL::abs(plattice[1][0]) <= 0.5) || (ANTL::abs(plattice[1][1]) <= 0.5)
         || ((plattice[1][0])*(plattice[1][1]) >= 0) ) {

           #ifdef DEBUG
           std::cout << "becomePrepared: Enter step 4,5,6" << std::endl;
           #endif
            //encompasses steps 3,4,5,6

      // STEP 3
          //computing a bound for the number of CF iterations

          //This picks the correct value for r, stores it in rR
          this->rR = Type(ceil(ANTL::sqrt( ANTL::abs(PType(ideal1.coeff_matrix[1][1])/PType(ideal1.coeff_matrix[2][2]))  )));

          if ( this->rR < ( (ideal1.coeff_matrix[0][0]) / (ideal1.coeff_matrix[2][2]) +1) ){
            this->rR =  (ideal1.coeff_matrix[0][0])/(ideal1.coeff_matrix[2][2]) +1;
          }
          //Note this is only valid if (a,b,c,d) is the reduced index form,
          //which  allows R = 2. The true bound which is needed is rR^2
          mul(this->rR, this->rR, 2);
          mul(this->rR, this->rR, this->rR); //squares rR




      // STEP 4
          NTL::clear(this->pm); NTL::set(this->pm1); NTL::set(this->qm); NTL::clear(this->qm1); m = Type(-2);
          div(this->alpha0, plattice[0][1] , plattice[0][0]);


      //STEP 5
          do {
            // At this point, the p
              ++m;
              this->a1 = Type(floor(this->alpha0));
              this->alpha1 = 1/(this->alpha0 - PType(this->a1));
              this->pNext = this->a1*this->pm1 + this->pm;  //(m+1)th
              this->qNext = this->a1*this->qm1 + this->qm; //(m+1)th

              //update variables
              this->pm = this->pm1;             // pm is now the mth term
              this->pm1 = this->pNext;          // pm1 is (m+1)th term

              this->qm = this->qm1;             // qm is the mth term
              this->qm1 = this->qNext;          // qm1 is the (m+1)th term
              this->alpha0 = this->alpha1;
              // std::cout <<  "m: " << m << "  p_m: " << pm <<  "  p_(m+1): " << pm1 <<std::endl;
            }
            while( !( (m%2 == 0) && (this->qm*this->qm1 > rR)  ) );

      //STEP 6. Transformation by [  pm1   -pm]
      //                          [ -qm1    qm]
          this->temp_pb[0][0] = PType(this->pm1)*plattice[0][0] - PType(this->qm1)*plattice[0][1];
          this->temp_pb[1][0] = PType(this->pm1)*plattice[1][0] - PType(this->qm1)*plattice[1][1];
          this->temp_pb[0][1] =  PType(this->qm)*plattice[0][1] -  PType(this->pm)*plattice[0][0];
          this->temp_pb[1][1] =  PType(this->qm)*plattice[1][1] -  PType(this->pm)*plattice[1][0];
          plattice[0][0] =this->temp_pb[0][0];
          plattice[1][0] =this->temp_pb[1][0];
          plattice[0][1] =this->temp_pb[0][1];
          plattice[1][1] =this->temp_pb[1][1];

          this->temp_lb[0][0] = this->pm1*ideal1.coeff_matrix[0][1] - this->qm1*ideal1.coeff_matrix[0][2];
          this->temp_lb[1][0] = this->pm1*ideal1.coeff_matrix[1][1] - this->qm1*ideal1.coeff_matrix[1][2];
          this->temp_lb[2][0] = this->pm1*ideal1.coeff_matrix[2][1] - this->qm1*ideal1.coeff_matrix[2][2];
          this->temp_lb[0][1] =  this->qm*ideal1.coeff_matrix[0][2] -  this->pm*ideal1.coeff_matrix[0][1];
          this->temp_lb[1][1] =  this->qm*ideal1.coeff_matrix[1][2] -  this->pm*ideal1.coeff_matrix[1][1];
          this->temp_lb[2][1] =  this->qm*ideal1.coeff_matrix[2][2] -  this->pm*ideal1.coeff_matrix[2][1];

          ideal1.coeff_matrix[0][1] = this->temp_lb[0][0];
          ideal1.coeff_matrix[1][1] = this->temp_lb[1][0];
          ideal1.coeff_matrix[2][1] = this->temp_lb[2][0];
          ideal1.coeff_matrix[0][2] = this->temp_lb[0][1];
          ideal1.coeff_matrix[1][2] = this->temp_lb[1][1];
          ideal1.coeff_matrix[2][2] = this->temp_lb[2][1];
      } //close the if of step 2

      //std::cout << "Prepared transformation 2" << std::endl;
      //std::cout <<  plattice[0][0] << "  " << plattice[0][1] <<  std::endl;
      //std::cout <<  plattice[1][0] << "  " << plattice[1][1] <<  std::endl;

      //STEP 7
      if ( !( ANTL::abs(plattice[1][0]) < 0.5 && ANTL::abs(plattice[1][1]) > 0.5 ) ){
          #ifdef DEBUG
          std::cout << "becomePrepared: Enter step 8,9,10" << std::endl;
          #endif
          //encompasses steps 8,9,10
          NTL::clear(this->pm); NTL::set(this->pm1); NTL::set(this->qm); NTL::clear(this->qm1); m = Type(-3);
          this->pminus = 0;
          this->alpha0 = - (plattice[1][1]) / plattice[1][0];
          //std::cout << "initial alpha " << alpha0<< std::endl;
          this->E = 2* ANTL::abs(plattice[1][0]);

          do {
              ++m;
              this->a1 = Type(floor(this->alpha0));        //a_(m+2)
              this->alpha1 = 1/(this->alpha0 - PType(this->a1));  //alpha_(m+3)
              this->pNext = this->a1*this->pm1 + this->pm;              //p_(m+2)
              this->qNext = this->a1*this->qm1 + this->qm;              //q_(m+2)


              //update variables
              this->pBefore = this->pminus;                 // This is p_(m-1)
              this->pminus = this->pm;                      // this is now p_m
              this->pm = this->pm1;                         // this is p_(m+1)
              this->pm1 = this->pNext;                      // p_(m+2)
              this->qBefore = this->qminus;                 // this is q_(m-1)
              this->qminus = this->qm;                      // q_m
              this->qm = this->qm1;                         // q_(m+1)
              this->qm1 = this->qNext;                      // q_(m+2)
              this->alpha0 = this->alpha1;                  // alpha_(m+3)
              //std::cout <<  "becomePrepared:: m: " << m << "  p_m: " << pminus <<  "  p_(m+1): " << pm << " alpha0 " << alpha0<<std::endl;
              //std::cout <<  "becomePrepared:: m: " << m << "  q_m: " << qminus <<  "  q_(m+1): " << qm <<std::endl;
          }while( !( PType(this->qm1) > E ) );
          if (m == -2){
            //std::cout << "ERROR: SCF ended on first round" << std::endl;
          }

          //std::cout << "CF step ended with q_{m+2} = " << qm1 << " which is larger than " << E << std::endl;
          //when this loop breaks, pminus = m, pm is m+1, and pm1 is m+2th position

          if( ANTL::abs(PType(pminus)*plattice[1][0] + PType(qminus)*plattice[1][1])>0.5 ){
            //std::cout << "Transform matrix" << std::endl;
            //std::cout <<  this->pm << "  " << this->pminus  <<  std::endl;
            //std::cout <<  this->qm << "  " << this->qminus <<  std::endl;

            #ifdef DEBUG
            std::cout << "becomePrepared: step 10: first transform" << std::endl;
            #endif
            this->temp_pb[0][0] = PType(this->pm)*plattice[0][0] + PType(this->qm)*plattice[0][1];
            this->temp_pb[1][0] = PType(this->pm)*plattice[1][0] + PType(this->qm)*plattice[1][1];
            this->temp_pb[0][1] = PType(this->pminus)*plattice[0][0] +  PType(this->qminus)*plattice[0][1];
            this->temp_pb[1][1] = PType(this->pminus)*plattice[1][0] +  PType(this->qminus)*plattice[1][1];

            plattice[0][0] = this->temp_pb[0][0];
            plattice[1][0] = this->temp_pb[1][0];
            plattice[0][1] = this->temp_pb[0][1];
            plattice[1][1] = this->temp_pb[1][1];


            this->temp_lb[0][0] = this->pm*ideal1.coeff_matrix[0][1] + this->qm*ideal1.coeff_matrix[0][2];
            this->temp_lb[1][0] = this->pm*ideal1.coeff_matrix[1][1] + this->qm*ideal1.coeff_matrix[1][2];
            this->temp_lb[2][0] = this->pm*ideal1.coeff_matrix[2][1] + this->qm*ideal1.coeff_matrix[2][2];
            this->temp_lb[0][1] =  this->qminus*ideal1.coeff_matrix[0][2] +  this->pminus*ideal1.coeff_matrix[0][1];
            this->temp_lb[1][1] =  this->qminus*ideal1.coeff_matrix[1][2] +  this->pminus*ideal1.coeff_matrix[1][1];
            this->temp_lb[2][1] =  this->qminus*ideal1.coeff_matrix[2][2] +  this->pminus*ideal1.coeff_matrix[2][1];

            ideal1.coeff_matrix[0][1] = this->temp_lb[0][0];
            ideal1.coeff_matrix[1][1] = this->temp_lb[1][0];
            ideal1.coeff_matrix[2][1] = this->temp_lb[2][0];
            ideal1.coeff_matrix[0][2] = this->temp_lb[0][1];
            ideal1.coeff_matrix[1][2] = this->temp_lb[1][1];
            ideal1.coeff_matrix[2][2] = this->temp_lb[2][1];

          }
          else{
            #ifdef DEBUG
            std::cout << "becomePrepared: Step10 second transform" << std::endl;
            #endif
            this->temp_pb[0][0] = PType(pminus)*plattice[0][0] + PType(qminus)*plattice[0][1];
            this->temp_pb[1][0] = PType(pminus)*plattice[1][0] + PType(qminus)*plattice[1][1];
            this->temp_pb[0][1] = PType(pBefore)*plattice[0][0] +  PType(qBefore)*plattice[0][1];
            this->temp_pb[1][1] = PType(pBefore)*plattice[1][0] +  PType(qBefore)*plattice[1][1];

            plattice[0][0] = this->temp_pb[0][0];
            plattice[1][0] = this->temp_pb[1][0];
            plattice[0][1] = this->temp_pb[0][1];
            plattice[1][1] = this->temp_pb[1][1];

            this->temp_lb[0][0] = this->pminus*ideal1.coeff_matrix[0][1] + this->qminus*ideal1.coeff_matrix[0][2];
            this->temp_lb[1][0] = this->pminus*ideal1.coeff_matrix[1][1] + this->qminus*ideal1.coeff_matrix[1][2];
            this->temp_lb[2][0] = this->pminus*ideal1.coeff_matrix[2][1] + this->qminus*ideal1.coeff_matrix[2][2];
            this->temp_lb[0][1] = this->qBefore*ideal1.coeff_matrix[0][2] + this->pBefore*ideal1.coeff_matrix[0][1];
            this->temp_lb[1][1] = this->qBefore*ideal1.coeff_matrix[1][2] + this->pBefore*ideal1.coeff_matrix[1][1];
            this->temp_lb[2][1] = this->qBefore*ideal1.coeff_matrix[2][2] + this->pBefore*ideal1.coeff_matrix[2][1];


            ideal1.coeff_matrix[0][1] = this->temp_lb[0][0];
            ideal1.coeff_matrix[1][1] = this->temp_lb[1][0];
            ideal1.coeff_matrix[2][1] = this->temp_lb[2][0];
            ideal1.coeff_matrix[0][2] = this->temp_lb[0][1];
            ideal1.coeff_matrix[1][2] = this->temp_lb[1][1];
            ideal1.coeff_matrix[2][2] = this->temp_lb[2][1];

          }

      }

      //std::cout << "Prepared transformation 3" << std::endl;
      //std::cout <<  plattice[0][0] << "  " << plattice[0][1] <<  std::endl;
      //std::cout <<  plattice[1][0] << "  " << plattice[1][1] <<  std::endl;
      //std::cout <<  ideal1.coeff_matrix[0][0] << "  " << ideal1.coeff_matrix[0][1] << "  " << ideal1.coeff_matrix[0][2] << std::endl;
      //std::cout <<  ideal1.coeff_matrix[1][0] << "  " << ideal1.coeff_matrix[1][1] << "  " << ideal1.coeff_matrix[1][2] << std::endl;
      //std::cout <<  ideal1.coeff_matrix[2][0] << "  " << ideal1.coeff_matrix[2][1] << "  " << ideal1.coeff_matrix[2][2] << std::endl;

      // Step 11
      this->dummy1 = Type(floor( plattice[0][1]/plattice[0][0] ));


      plattice[0][1] -= PType(this->dummy1)*plattice[0][0];
      plattice[1][1] -= PType(this->dummy1)*plattice[1][0];

      ideal1.coeff_matrix[0][2] -= this->dummy1*ideal1.coeff_matrix[0][1];
      ideal1.coeff_matrix[1][2] -= this->dummy1*ideal1.coeff_matrix[1][1];
      ideal1.coeff_matrix[2][2] -= this->dummy1*ideal1.coeff_matrix[2][1];

      ideal1.normalize();


}

//template <typename Type, typename PType>
//void VoronoiMethods<Type, PType>::voronoi_step(CubicIdeal<Type, PType> & ideal1){
//
//}

#endif // guard
