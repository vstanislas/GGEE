#include <RcppArmadillo.h>
#include <vector>
#include <iostream>

// [[Rcpp::plugins(cpp11)]]


//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//[[Rcpp::export]]
List  SVD(arma::mat x){ 
    arma::mat U ;
    arma::vec s ;
    arma::mat V ;
	
	svd(U, s, V, x) ;
	
	List ret ;
	ret["U"] = U ;
	ret["s"] = s ;
	ret["V"] = V ;
  	
  	return(ret) ;
 }
 


 //[[Rcpp::export]]
NumericMatrix cbind1 (NumericVector x, NumericVector y){
  NumericMatrix out (x.size(), 2);
  out(_,0) = x; 
  out(_,1)=y;
  return out;
}

 
 //[[Rcpp::export]]
NumericMatrix cbind2 (NumericMatrix a, NumericVector b) {
  int acoln = a.ncol();
  int bcoln = 1;
  NumericMatrix out = no_init_matrix(a.nrow(), acoln + bcoln);
  for (int j = 0; j < acoln + bcoln; j++) {
    if (j < acoln) {
      out(_, j) = a(_, j);
    } else {
      out(_, j) = b ;
    }
  }
  return out;
}


//[[Rcpp::export]]
List  IntProd(NumericMatrix g1, NumericMatrix g2, int g1num, int g2num){ 

	int p1 = g1.ncol() ;
	int p2 = g2.ncol() ;
  	
    NumericMatrix W(g1.nrow(),0) ;
	CharacterVector names;
	
  	for(int i = 0; i < p1; i++){
    	for(int j = 0; j < p2; j++){
      		
      		NumericVector snp1 = g1(_, i) ;
      		NumericVector snp2 = g2(_, j) ;
      		int n = snp1.size() ;
      		NumericVector prod(n) ;
      	
      		for(int k =0; k<n; k++){
      			int prod_temp = snp1[k] * snp2[k] ;
				prod[k] = prod_temp ;
     		}
     		
      		// Rcout << prod<<"\n" ;
      		W = cbind2(W, prod) ;
      		std::string result = "pair." + std::to_string(g1num) + "." + std::to_string(g2num) + ".snp." + std::to_string(i+1) + "."  + std::to_string(j+1) ;
      		names.push_back(result) ;
      	}
  	}
  
  	List ret ;
  	ret["W"]=W ;
  	ret["names"] = names ;
  	
	return ret;
}


            
          
