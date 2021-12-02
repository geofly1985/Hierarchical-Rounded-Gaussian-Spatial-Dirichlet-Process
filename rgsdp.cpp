//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <iostream>
#include <algorithm> // unique
#include </Users/lix233/Desktop/JRSS-C/code/dp.h> //wsample

using namespace std;
using namespace arma;
using namespace Rcpp;



//[[Rcpp::export]]
double arma_sub(rowvec x, uvec pos) {

    rowvec z = pow(x,2);
    double s = sum(z.elem(pos));  // sum Xj^2

    return s;
}



//[[Rcpp::export]]
vec mvrnorm(vec M, mat S) {
  int n = M.n_rows;
  mat e = randn(n);
  return vectorise(M + chol(S).t()*e);
}



//[[Rcpp::export]]
mat uniqueRows(mat Y) { // This works and is fast, but not certain of accuracy
  int n = Y.n_rows;
  int k = Y.n_cols;
  vec x = Y * randu(k);
  vec idx = zeros<vec>(n);

  map<double, int> uRow;
  for (int i=0; i<n; i++) {
    if ( uRow.find(x[i]) == uRow.end() ) { // key doesn't exist
      uRow[x[i]] = 1;
      idx[i] = 1;
    } else { // key found
      uRow[x[i]] = uRow[x[i]] + 1;
    }
  }

  return Y.rows( find(idx == 1) );
}



//[[Rcpp::export]]
uvec matchRows (mat X, vec v) {
  uvec uv;
  vec vv = zeros<vec>(X.n_rows);
  mat out;

  for (int i=0; i<X.n_rows; i++) {
    vv[i] = 1 * all(vectorise(X.row(i)) == v);
  }

  uv = find(vv == 1);
  return uv;
}



//[[Rcpp::export]]
double ldmvnorm(vec x, vec m, mat S) {
  int n = x.size();
  double ldet, val, sign;
  vec out;
  vec b = x - m;
  log_det(val,sign,S);
  ldet = val;
  out = -n/2*log(2*pi) -.5*ldet -.5*b.t()*S.i()*b;

  return out[0];
}


//[[Rcpp::export]]
arma::cube cube_trans(NumericVector B_) {
  IntegerVector dimB = B_.attr("dim");
  arma::cube D(B_.begin(), dimB[0], dimB[1], dimB[2]);

  return D;

}


//[[Rcpp::export]]
mat update_theta(vec X, mat Y, mat D, mat W, double rho, mat Theta, double sig2, double tau2, double v, vec Beta, NumericVector RR){

  cube R = cube_trans(RR);
  double ldet, val, sign;
  double ldetH, val2, sign2;

  mat Theta_new = Theta;
  mat H = (D - rho*W).i();
  mat H_inv = D - rho*W;
  int T = Y.n_rows;
  int n = Y.n_cols;
  mat Theta_new_c = zeros(T,n);

  mat Lambda, Lambda_j;
  vec mu, mu_j, mean_all;
  mat MU, sum_j;
  double logq0, logqj, Tj, ind;
  mat Theta_xt, Theta_xt_star, Theta_star;
  vec Tind = linspace(0, T-1, T);
  vec lprobs, probs;
  uvec vv, inds;
  int J, J_xt, Tj_xt;
  mat In = eye<mat>(n,n);

  mat tmp, tmp1;


  log_det(val2, sign2, H);
  ldetH = val2;




  for (int t=0; t<T; t++) {

    Lambda = ((pow(X[t],2)/tau2)*In + (1/sig2)*H_inv).i();
    log_det(val, sign, Lambda);
    ldet = val;


    mu = vectorise(Y.row(t)) - R.slice(t).t()*Beta;

    tmp = log(v) + 0.5*ldet - 0.5*ldetH - 0.5*n*log(2*pi*tau2*sig2)- 0.5/tau2*(mu.t()*(In- pow(X[t],2)/tau2*Lambda)*mu);
    logq0 = tmp(0,0);

    Theta_xt = Theta.rows(find(Tind != t));

    Theta_xt_star = uniqueRows(Theta_xt);

    J_xt = Theta_xt_star.n_rows;

    lprobs = zeros(J_xt+1);
    probs = zeros(J_xt+1);

    lprobs[J_xt] = logq0;


      for (int j=0; j<J_xt; j++) {

        vv = matchRows(Theta_xt,vectorise( Theta_xt_star.row(j) ));

        Tj_xt = vv.size();

        tmp1 = -n/2*log(2*pi) -n/2*log(tau2)- 0.5/tau2*( (mu-X[t]*(vectorise(Theta_xt_star.row(j)))).t()*(mu-X[t]*(vectorise(Theta_xt_star.row(j)))));
        logqj = tmp1(0, 0);
        lprobs[j] = log(Tj_xt) + logqj;

      }


    probs = exp(lprobs - max(lprobs));
    ind = wsample(linspace(0,J_xt,J_xt+1), probs);
    if (ind == J_xt) {
      Theta_new.row(t) = reshape(mvrnorm((X[t]/tau2)*Lambda*mu,Lambda), 1, n);
    } else {
      Theta_new.row(t) = Theta_xt_star.row(ind);
    }


  }


// update Theta after leave one out step
  Theta_star = uniqueRows(Theta_new);
  J = Theta_star.n_rows;

  for (int j=0;j<J;j++) {

    inds = matchRows(Theta_new, vectorise(Theta_star.row(j)));

    Tj = arma_sub(X.t(),inds);

    // jth cluster Covariance matrix
    Lambda_j = (Tj/tau2*In + H_inv/sig2).i();

    // calculate sum_j
    sum_j = zeros(n,1);
    for (int k=0; k<inds.n_rows; k++){
    sum_j = sum_j + X[inds[k]]*(vectorise(Y.row(inds[k])) - R.slice(inds[k]).t()*Beta);
    }

    // jth cluster mean vector
    mu_j = Lambda_j/tau2*sum_j;


    Theta_new.each_row(inds) = reshape(mvrnorm(mu_j, Lambda_j), 1, n);


  }



  mean_all = sum(Theta_new, 1)/n;

  for (int k=0;k<T;k++){
    Theta_new_c.row(k) = Theta_new.row(k) - ones(1,n)*(mean_all.row(k)[0]);
  }




  return Theta_new_c;

}
