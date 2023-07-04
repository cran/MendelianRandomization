#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
using namespace Rcpp;
using namespace arma;

double l_i_c(const arma::vec& b_exp_i,
           const arma::vec& b_out_i,
           const arma::mat& Sig_inv_i,
           const arma::vec& b_t_i,
           const arma::vec& theta_t,
           const arma::vec& r_vec_t_i)
{
  arma::vec beta = arma::join_vert(b_exp_i, b_out_i);
  arma::vec bvec = arma::join_vert(b_t_i, dot(theta_t, b_t_i) + r_vec_t_i);
  double l =  - 1.0/2.0 * arma::as_scalar((beta-bvec).t() * Sig_inv_i * (beta-bvec));
  return l;
}

double loglik_c(const arma::mat& b_exp, const arma::mat& b_out,
                const Rcpp::List& Sig_inv_l,
                const arma::mat& b_t,
                const arma::vec& theta_t,
                const arma::vec& r_vec_t){
  int p = b_exp.n_rows;
  arma::vec l_i_vec(p);
  for(int i=0; i < p; ++i){
   l_i_vec(i)  = l_i_c(b_exp.row(i).t(),b_out.row(i), Sig_inv_l[i],b_t.row(i).t(),theta_t,r_vec_t.row(i));
  }
  double l = sum(l_i_vec);
  return l;
}

Rcpp::List WK_func(const Rcpp::List& Sig_inv_l)
  {
  const int p = Sig_inv_l.length(), L = as<arma::mat>(Sig_inv_l[0]).n_rows;
  arma::mat WK1mat(p,L);
  arma::vec WK1vec(p);
  for(int i = 0; i < p; i++){
    WK1mat.row(i) = as<arma::mat>(Sig_inv_l[i]).row(L-1);
    WK1vec(i) = as<arma::mat>(Sig_inv_l[i])(L-1,L-1);
  }
  return Rcpp::List::create(Rcpp::Named("WK1vec") = WK1vec,
                            Rcpp::Named("WK1mat") = WK1mat);
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

Rcpp::List MVcML_estimate_c(const arma::mat& b_exp, const arma::mat& b_out,
                          const Rcpp::List& Sig_inv_l,
                          const arma::vec& WK1vec,
                          const arma::mat& WK1mat,
                          const unsigned int K,
                          const arma::vec& initial_theta,
                          const arma::mat& initial_mu,
                          const unsigned int maxit = 100,
                          const double thres = 1e-4)
{

  // initialize
  const int p = b_exp.n_rows, k = b_exp.n_cols;
  arma::vec theta = initial_theta, theta_old = theta -1;
  arma::mat mu_vec = initial_mu;
  arma::vec v_importance(p), v_bg(p), beta_star_i, B;
  arma::mat r_t(p,1);
  arma::mat A, temp, beta_star;
  unsigned int ite_ind;
  //
  for (ite_ind = 0; ite_ind < maxit && sum(abs(theta - theta_old)) > thres; ++ite_ind){
    theta_old = theta;

    if (K > 0) {
      A = arma::eye(p,p);
      A.diag() = 1 / WK1vec;
      beta_star = arma::join_horiz(b_exp-mu_vec,b_out-(mu_vec * theta));
      temp = WK1mat * beta_star.t();
      B = temp.diag();
      r_t = A * B;
      for(int i = 0; i < p; ++i){
        double l_r = l_i_c(b_exp.row(i).t(), b_out.row(i), Sig_inv_l[i], mu_vec.row(i).t(), theta, r_t.row(i));
        double l_r0 = l_i_c(b_exp.row(i).t(), b_out.row(i), Sig_inv_l[i], mu_vec.row(i).t(), theta, zeros(1,1));
        v_importance(i) = l_r-l_r0;
      }

      arma::uvec nonzero_bg_ind = sort_index(v_importance, "descend");
      v_bg = arma::zeros(p);
      v_bg(nonzero_bg_ind.head(K)) = r_t(nonzero_bg_ind.head(K));
    }
    else{
      v_bg = arma::zeros(p);
    }
    arma::mat theta_cp = theta * theta.t();
    for(int i = 0; i < p; ++i){
          beta_star_i = arma::join_vert(b_exp.row(i).t(), b_out.row(i) - v_bg.row(i));
          arma::mat W_i = Sig_inv_l[i];
          B = W_i.rows(0, k-1) * beta_star_i +
            arma::as_scalar(W_i.row(k) * beta_star_i) * theta;
          A = W_i(span(0,k-1), span(0,k-1)) +  arma::as_scalar(W_i(k,k)) * theta_cp +
            W_i(span(0,k-1),k) * theta.t() + theta * W_i(k,span(0,k-1));

          mu_vec.row(i) = solve(A,B).t();
        }

    A = zeros(k,k);
    B = zeros(k);
    for(int i = 0; i < p; ++i){
      arma::mat W_i = Sig_inv_l[i];
      A = A + arma::as_scalar(W_i(k,k)) * mu_vec.row(i).t() * mu_vec.row(i);
      beta_star_i = arma::join_vert(b_exp.row(i).t()-mu_vec.row(i).t(),b_out.row(i)-v_bg.row(i));
      B = B + arma::as_scalar(W_i.row(k) * beta_star_i) * mu_vec.row(i).t();
    }

    theta = solve(A,B);

  }
  bool Conv = (ite_ind >= maxit);
  double Neg_l = -loglik_c(b_exp,b_out,Sig_inv_l,mu_vec,
                          theta,v_bg);

  return Rcpp::List::create(Rcpp::Named("theta") = theta,
                           Rcpp::Named("b_vec") = mu_vec,
                           Rcpp::Named("r_vec") = v_bg,
                           Rcpp::Named("l") = Neg_l,
                           Rcpp::Named("Converge") = Conv);
}

Rcpp::List MVcML_estimate_random_c(const arma::mat& b_exp, const arma::mat& b_out,
                                   const arma::mat& se_bx,
                                   const Rcpp::List& Sig_inv_l,
                                   const arma::vec& WK1vec,
                                   const arma::mat& WK1mat,
                                   const unsigned int K,
                                   const unsigned int n,
                                   const unsigned int random_start = 1,
                                   const double min_theta_range = -0.5,
                                   const double max_theta_range = 0.5,
                                   const unsigned int maxit = 100,
                                   const double thres = 1e-4)
{
// initialize
  const int p = b_exp.n_rows, k = b_exp.n_cols;
  arma::mat theta_v_RandomCandidate(k,random_start), invalid_RandomCandidate(p,random_start);
  arma::vec l_v_RandomCandidate(random_start), Conv_v_RandomCandidate(random_start);
  arma::vec initial_theta(k);
  arma::mat initial_mu(p,k);
  for(unsigned int random_ind = 0; random_ind < random_start; ++random_ind)
  {
    if(random_ind == 0)
    {
      initial_theta = zeros(k);
      // initial_mu = zeros(p,k);
      initial_mu = b_exp;
    } else {
      initial_theta = runif(k,  min_theta_range,  max_theta_range); // TO-DO: option to change initial value for EACH theta
      arma::mat b_rand;
      b_rand = Rcpp::rnorm(p*k,0.0,1.0);
      b_rand.reshape(size(b_exp));
      initial_mu = b_exp + b_rand % se_bx;
      // initial_mu = zeros(p,k);
    }

    Rcpp::List MLE_result = MVcML_estimate_c(b_exp,b_out,
                              Sig_inv_l,WK1vec,WK1mat,
                              K,initial_theta,
                              initial_mu,
                              maxit,
                              thres);

    double Neg_l = -loglik_c(b_exp,b_out,Sig_inv_l,MLE_result["b_vec"],
                               MLE_result["theta"],MLE_result["r_vec"]);

    theta_v_RandomCandidate.col(random_ind) = as<arma::vec>(MLE_result["theta"]);
    l_v_RandomCandidate(random_ind) = Neg_l;
    Conv_v_RandomCandidate(random_ind) = as<int>(MLE_result["Converge"]);
    invalid_RandomCandidate.col(random_ind) = as<arma::vec>(MLE_result["r_vec"]);

  }
  uvec q1 = find(Conv_v_RandomCandidate == 0);
  if (q1.n_elem == 0){
    return Rcpp::List::create(Rcpp::Named("Converge") = 1);
  }
  int min_neg_l = l_v_RandomCandidate(q1).eval().index_min();

  arma::vec theta_est = theta_v_RandomCandidate.cols(q1).eval().col(min_neg_l);
  double l_est = (l_v_RandomCandidate(q1).eval())(min_neg_l);
  arma::vec r_est = invalid_RandomCandidate.cols(q1).eval().col(min_neg_l);
  return Rcpp::List::create(Rcpp::Named("theta") = theta_est,
                            Rcpp::Named("l") = l_est,
                            Rcpp::Named("r_vec") = r_est,
                            Rcpp::Named("Converge") = 0);
}

//' MVMRcML method with Data Perturbation
//'
//' This is the internal MVMRcML-BIC function of mr_mvcML.
//'
//' @param b_exp A m*L matrix of SNP effects on the exposure variable.
//' @param b_out A m*1 matrix of SNP effects on the outcome variable.
//' @param se_bx A m*L matrix of standard errors of \code{b_exp}.
//' @param Sig_inv_l A list of the inverse of m covariance matrices.
//' @param n The smallest sample size of the L+1 GWAS dataset.
//' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
//' @param random_start Number of random start points, default is 1.
//' @param min_theta_range The lower bound of the uniform distribution for each initial value for theta generated from.
//' @param max_theta_range The upper bound of the uniform distribution for each initial value for theta generated from.
//' @param maxit Maximum number of iterations for each optimization, default is 100.
//' @param thres Threshold for convergence criterion.
//'
//' @import Rcpp
//' @import RcppArmadillo
//' @return A list
//' \describe{
//' \item{BIC_theta}{Estimated causal effect from MVMR-cML-BIC}
//' \item{BIC_invalid}{Invalid IVs selected by MVMR-cML-BIC}
//' \item{l_vec}{A vector of negative log-likelihood corresponding to each \code{K}.}
//' \item{K_vec}{A vector of candidate K's}
//' \item{theta_vec}{A matrix of causal parameter estimates, each column corresponds to a candidate \code{K}.}
//' \item{Conv_vec}{A vector of successful convergence indicators corresponding to each \code{K}.}
//' \item{Converge}{Indicator of successful convergence, 0 means success, 1 means failure.}
//' \item{BIC_vec}{Data perturbation with successful convergence}
//' \item{Khat}{The length of \code{BIC_invalid}.}
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List MVmr_cML(const arma::mat& b_exp, const arma::mat& b_out,
                      const arma::mat& se_bx,
                      const Rcpp::List& Sig_inv_l,
                      const unsigned int n,
                      Rcpp::NumericVector K_vec = Rcpp::NumericVector::create(),
                      const unsigned int random_start = 1,
                      const double min_theta_range = -0.5,
                      const double max_theta_range = 0.5,
                      const unsigned int maxit = 100,
                      const double thres = 1e-4)
{
  //initialize
  const int p = b_exp.n_rows, L = b_exp.n_cols;
  arma::vec K_vec1;
  if(K_vec.size() == 0){
    K_vec1 = linspace(0, p - L - 1, p - L);
  }else{
    K_vec1 = as<arma::vec>(K_vec);
  }
  arma::vec WK1vec = as<arma::vec>(WK_func(Sig_inv_l)["WK1vec"]);
  arma::mat WK1mat = as<arma::mat>(WK_func(Sig_inv_l)["WK1mat"]);

  const int K_vec_len = K_vec1.n_elem;
  arma::mat rand_theta(L,K_vec_len), invalid_mat(p,K_vec_len);
  arma::vec rand_l(K_vec_len), Conv_l(K_vec_len), BIC_vec;
  for(int i=0; i < K_vec_len; ++i)
  {
    int K_value = K_vec1(i);
    Rcpp::List rand_res = MVcML_estimate_random_c(b_exp,b_out,se_bx,Sig_inv_l,
                                     WK1vec, WK1mat,
                                     K_value, n,
                                     random_start, min_theta_range, max_theta_range, maxit,thres);
    if(as<int>(rand_res["Converge"]) == 1){
      Conv_l(i) = 1;
      continue;
    }
    rand_theta.col(i) = as<arma::vec>(rand_res["theta"]);
    rand_l(i)= rand_res["l"];
    Conv_l(i) = 0;
    invalid_mat.col(i) = as<arma::vec>(rand_res["r_vec"]);
  }
  uvec q1 = find(Conv_l == 0);
  if(q1.n_elem==0){
    return Rcpp::List::create(Rcpp::Named("Converge") = 1);
  }
  BIC_vec = log(n) * K_vec1(q1) + 2 * rand_l(q1);

// cML-BIC
  unsigned int min_ind = BIC_vec.index_min();
  uvec BIC_invalid = arma::find(invalid_mat.cols(q1).eval().col(min_ind));

  arma::vec BIC_theta = rand_theta.cols(q1).eval().col(min_ind);

      Rcpp::List res;
      res["BIC_theta"] = BIC_theta;
      res["BIC_invalid"] = BIC_invalid + 1;
      res["l_vec"] = rand_l;
      res["K_vec"] = K_vec1;
      res["theta_vec"] = rand_theta;
      res["Conv_vec"] = Conv_l;
      res["Converge"] = 0;
      res["BIC_vec"] = BIC_vec;
      res["Khat"] = BIC_invalid.n_elem;
      return res;

}

//' MVMRcML method with Data Perturbation
//'
//' This is the internal MVMRcML-DP function of mr_mvcML.
//'
//' @param b_exp A m*L matrix of SNP effects on the exposure variable.
//' @param b_out A m*1 matrix of SNP effects on the outcome variable.
//' @param se_bx A m*L matrix of standard errors of \code{b_exp}.
//' @param Sig_inv_l A list of the inverse of m covariance matrices.
//' @param n The smallest sample size of the L+1 GWAS dataset.
//' @param K_vec Sets of candidate K's, the constraint parameter representing number of invalid IVs.
//' @param random_start Number of random start points, default is 1.
//' @param num_pert Number of perturbation, default is 100.
//' @param min_theta_range The lower bound of the uniform distribution for each initial value for theta generated from.
//' @param max_theta_range The upper bound of the uniform distribution for each initial value for theta generated from.
//' @param maxit Maximum number of iterations for each optimization, default is 100.
//' @param thres Threshold for convergence criterion.
//'
//' @import Rcpp
//' @import RcppArmadillo
//' @return A list
//' \describe{
//' \item{BIC_theta}{Estimated causal effect from MVMR-cML-BIC}
//' \item{BIC_invalid}{Invalid IVs selected by MVMR-cML-BIC}
//' \item{BIC_DP_theta}{Estimated causal effect from MVMR-cML-DP }
//' \item{BIC_DP_se}{Estimate standard error for \code{BIC_DP_theta}}
//' \item{eff_DP_B}{Data perturbation with successful convergence}
//' \item{DP_ninvalid}{A vector of the number of selected invalid IVs by MVMRcML-BIC in each data perturbation.}
//' }
//'
//' @export
// [[Rcpp::export]]
Rcpp::List MVmr_cML_DP(const arma::mat& b_exp, const arma::mat& b_out,
                      const arma::mat& se_bx,
                      const Rcpp::List& Sig_inv_l,
                      const unsigned int n,
                      Rcpp::NumericVector K_vec = Rcpp::NumericVector::create(),
                      const unsigned int random_start = 1,
                      const unsigned int num_pert = 100,
                      const double min_theta_range = -0.5,
                      const double max_theta_range = 0.5,
                      const unsigned int maxit = 100,
                      const double thres = 1e-4)
{
  //initialize
  const int p = b_exp.n_rows, L = b_exp.n_cols;
  //arma::vec K_vec1;
  arma::mat theta_v(L,num_pert);
  Rcpp::List cML_res = MVmr_cML(b_exp, b_out, se_bx, Sig_inv_l, n, K_vec, random_start, min_theta_range, max_theta_range, maxit, thres);
  arma::mat b_exp_new(size(b_exp)), b_out_new(size(b_out));
  arma::vec theta_b, Conv_v(num_pert), K_v(num_pert);
  for(unsigned int pt_ind = 0; pt_ind < num_pert; ++pt_ind)
    {
    for(int i = 0; i < p; ++i)
      {
        arma::mat sig = arma::inv(as<arma::mat>(Sig_inv_l[i]));
        arma::mat epis = mvrnormArma(1, zeros(L+1), sig);
        b_exp_new.row(i) = b_exp.row(i) + epis.head_cols(L);
        b_out_new.row(i) = b_out.row(i) + epis.col(L);
      }
      Rcpp::List cML_res_b = MVmr_cML(b_exp_new, b_out_new, se_bx, Sig_inv_l, n, K_vec, random_start, min_theta_range, max_theta_range, maxit, thres);
      if(as<int>(cML_res_b["Converge"])==1){
        Conv_v(pt_ind) = 1;
        continue;
      }
      theta_b = as<arma::vec>(cML_res_b["BIC_theta"]);
      theta_v.col(pt_ind) = theta_b;
      Conv_v(pt_ind) = 0;
      K_v(pt_ind) = as<int>(cML_res_b["Khat"]);
    }
  uvec q1 = find(Conv_v == 0);
  // cML-BIC-DP
  arma::vec BIC_DP_theta = mean(theta_v.cols(q1).eval(), 1);
  arma::vec BIC_DP_se = stddev(theta_v.cols(q1).eval(), 0, 1);

  Rcpp::List res;
  res["BIC_theta"] = cML_res["BIC_theta"];
  res["BIC_invalid"] = cML_res["BIC_invalid"];
  res["BIC_DP_theta"] = BIC_DP_theta;
  res["BIC_DP_se"] = BIC_DP_se;
  res["eff_DP_B"] = q1.n_elem;
  res["DP_ninvalid"] = K_v;
  return res;

}




