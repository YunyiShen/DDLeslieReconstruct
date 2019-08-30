// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h> // to use sparse matrix
#include <climits>
using namespace Rcpp;

///get Leslie matrix from survival etc.
// [[Rcpp::export]]
arma::mat getLeslie(const arma::mat& Surv, const arma::mat& Fec, double SRB, bool minus1){
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
	if (minus1){
		Tr.diag() = Tr.diag()-1;
	}
	return(Tr);
}

///Calculate the density dependency
// [[Rcpp::export]]
arma::mat DD(bool global, arma::mat& Xn, arma::mat& E0, double aK0, const arma::mat& midP, bool null,const IntegerVector& nage){
  E0 = E0/sum(E0);
  arma::mat D;
  if(global){
    arma::mat den = ( 1 + (aK0) * (sum(Xn)-midP) );
    D = (1-null)*den + null;
  }
  else{
    D = 1-(1-null)* aK0 * (Xn-midP);
  }
  return(D);
}
