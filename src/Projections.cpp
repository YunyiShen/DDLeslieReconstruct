// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <climits>
using namespace Rcpp;

///get Leslie matrix from survival etc.
// [[Rcpp::export]]
arma::mat getLeslieCpp(const arma::mat& Surv, const arma::mat& Fec, const double& SRB){
	int nage_female = Fec.n_elem;
	int nage_male = Surv.n_elem - nage_female;
	arma::mat Tr(nage_female + nage_male , nage_male + nage_female);
	Tr.zeros();
	Tr.row(0).cols(0,nage_female-1) = SRB * Fec.t();
	Tr.row(nage_female).cols(0,nage_female-1) = (1 - SRB) * Fec.t();
	Tr.submat(1,0,nage_female-1,nage_female-2).diag() = Surv.rows(0,nage_female-2).t();
	Tr(nage_female-1,nage_female-1) = Surv(nage_female-1,0);
	Tr.submat(nage_female + 1,nage_female,nage_male + nage_female-1,nage_male + nage_female-2).diag() = Surv.rows(nage_female,nage_male + nage_female-2).t();
	Tr(nage_male+nage_female-1,nage_male+nage_female-1) = Surv(nage_male+nage_female-1,0);
	return(Tr);
}

///Calculate the density dependency
//[[Rcpp::export]]
arma::mat DD(const bool& global, const arma::mat& Xn,const arma::mat & aK0, const arma::mat& midP, const bool& null){
  //E0 = E0/sum(E0);// This was done in main projector
  arma::mat D;
  if(global){
    arma::mat den = ( 1 + (aK0) * (sum(Xn)-midP) );
    D = (1-null)*den + null;
  }
  else{
    D = 1-(1-null)* aK0 % (Xn-midP);
  }
  return(D);
}

///Helper function for a single year projection, inner function, export for test.
//[[Rcpp::export]]
arma::mat ProjectHarvest_helperCpp(const arma::mat& data_n,const arma::mat& Surv, const arma::mat& Fec,const double& SRB,const arma::mat& H_n, const arma::mat& H_np1,bool global, const List& aK0,const bool & null){
	arma::mat X_n1 = (1-H_n) % (data_n/H_n);
	arma::mat D_bir = DD(global, X_n1, aK0[0], aK0[2] ,null);
	arma::mat D_dea = DD(global, X_n1, aK0[1], aK0[2] ,null);
	return(H_np1 % (getLeslieCpp(Surv % D_dea, Fec % D_bir, SRB)*X_n1));
	
}

///main projection function
//[[Rcpp::export]]
arma::mat ProjectHarvestCpp(const arma::mat& Surv,const arma::mat& Harvpar,const arma::mat& Fec, const arma::mat& SRB, const List& aK0, const bool& global, const bool& null, const arma::mat& bl ,const int& period, const IntegerVector& nage){
	arma::mat Harvest(sum(nage),period+1);
	Harvest.col(0) = bl;
	//E0 = E0/(sum(E0));// need to check whether there is one in R call rather than here.
	for(int i = 1; i<period + 1; i++){
		Harvest.col(i) = ProjectHarvest_helperCpp(Harvest.col(i-1),Surv.col(i-1),Fec.col(i-1),(SRB(0,i-1)), Harvpar.col(i-1),Harvpar.col(i),global, aK0,null);
	}
	return(Harvest);
}

///get Aerial count
//[[Rcpp::export]]
arma::mat getAerialCountCpp(const arma::mat& Harv, const arma::mat& H, const arma::mat& A){
  return((sum((1/H-1) % Harv))%A);
}

///Misc



