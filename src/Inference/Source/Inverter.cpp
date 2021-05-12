#include "../Include/Inverter.h"


void Inverse_Base::print_for_debug(void) const {
  
  Rprintf( "Inverse computed: %d \n ", inverse_computed);
  if(inverse_computed){
    Rprintf( "Inverse E_inv (only some samples): \n");
    for (UInt i=0; i<10; i++){
      Rprintf( "E_inv( %d, %d):  %f \n", 10*i, 20*i, E_inv(10*i,20*i));
    } 
  }
  return;
};

void Inverse_Exact::Compute_Inv(){
  if(!this->inverse_computed){
    E_inv=this->E_decp->solve(MatrixXr::Identity(this->Ep->rows(),this->Ep->cols())); //Solve directly the system for an identity matrix
    this->inverse_computed=true;
  }
  
  return;
};

void Inverse_Non_Exact::Compute_Inv(){
  if(!this->inverse_computed){
    
    Eigen::BiCGSTAB<SpMat> Iterative_Solver;
    Iterative_Solver.compute(*(this->Ep));
    E_inv=Iterative_Solver.solve(MatrixXr::Identity(this->Ep->rows(),this->Ep->cols())); //Solve directly the system for an identity matrix (Iterative solution)
    
    this->inverse_computed=true;
  }
  
  return;
};
