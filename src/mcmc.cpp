// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// ====================================Supporting functions====================================
Rcpp::NumericVector grad(Function Model, List Data, NumericVector par, double h) {
  NumericMatrix mat(par.length(), par.length());
  NumericVector out(par.length());
  double        f_x = as<double>(as<List>(Model(par, Data))["LP"]);
  
  for (int i = 0; i < mat.ncol(); i++) {
    mat(i, _) = par;
    mat(i, i) = mat(i, i) + h;
    out[i] = (as<double>(as<List>(Model(mat(i, _), Data))["LP"]) - f_x) / h;
  }
  
  return out;
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  arma::mat m = Y * arma::chol(sigma);
  return m;
}

Rcpp::NumericVector HARproposal(NumericVector par) {
  
  NumericVector mu(par.length());
  arma::mat     Sigma = Rcpp::as<arma::mat>(NumericMatrix::diag(par.length(), 1));
  Rcpp::NumericVector theta = as<NumericVector>(wrap(mvrnormArma(1, mu, Sigma)));
  Rcpp::NumericVector d = theta / sqrt(sum(theta * theta));
  double              u = as<double>(runif(1));
  Rcpp::NumericVector prop = par + (u * d);
  
  return prop;
}

Rcpp::NumericVector SHARproposal(NumericVector par, double h, List Data, Function Model, NumericVector gr) {
  
  NumericVector mu(par.length());
  arma::mat     Sigma = Rcpp::as<arma::mat>(NumericMatrix::diag(par.length(), 1));
  Rcpp::NumericVector theta = as<NumericVector>(wrap(mvrnormArma(1, mu, Sigma)));
  Rcpp::NumericVector d = theta / sqrt(sum(theta * theta));
  Rcpp::NumericVector stdgr = gr / sqrt(sum(gr * gr));
  double              u = as<double>(runif(1));
  Rcpp::NumericVector prop = par + (u * (stdgr + d));

  return prop;
}

// ====================================MCMC Algorithms====================================

// [[Rcpp::export]]
SEXP harmwg(Function Model, List Data, int Iterations, int Status,
            int Thinning, double Acceptance, NumericMatrix Dev,
            int LIV, NumericMatrix Mon, List Mo0, NumericMatrix thinned) {
  
  // Initial settings
  int t_iter = 0, fins = 0, mcols = Mon.ncol();
  double alpha = 0;
  Rcpp::List Mo1 = clone(Mo0);
  RNGScope scope;
  
  // Run MCMC algorithm
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 <<
        ",   Proposal: Componentwise,   LP: " <<
          floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
      t_iter = floor((iter) / Thinning) + 1;
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
    // Random-Scan Componentwise Estimation
    Rcpp::NumericVector u = runif(LIV),      // Acceptance threshold for each parameter
                        z = rnorm(LIV);      // New Proposed value
    Rcpp::IntegerVector LIVseq = Rcpp::Range(0, LIV - 1), // Indexes to sample
      s = Rcpp::RcppArmadillo::sample(LIVseq, LIV, // Sample of indexes
                                      false, NumericVector::create());
    // Propose and evaluate new values per parameter
    for (int j = 0; j < LIV; j++) {
      // Propose new parameter values
      Rcpp::List Mo0_ = clone(Mo0);
      Rcpp::NumericVector prop = Mo0_["parm"];
      Rcpp::NumericVector prop1 = HARproposal(prop);
      prop[s[j]] = prop1[s[j]];
      // Log-Posterior of the proposed state
      Rcpp::List Mo1 = Model(prop, Data);
      fins = ::R_finite(Mo1["LP"]) + ::R_finite(Mo1["Dev"]);
      for (int m = 0; m < mcols; m++) {
        fins += ::R_finite(as<Rcpp::NumericVector>(Mo1["Monitor"])[m]);
      }
      if (fins < (mcols + 2)) Mo1 = Mo0;
      // Accept/Reject
      double LP0 = Mo0_["LP"];
      double LP1 = Mo1["LP"];
      alpha = exp(LP1 - LP0);
      if (u[s[j]] < alpha) {
        Mo0 = Mo1;
        Acceptance += 1.0 / LIV;
      }
    }
    if ((iter + 1) % Thinning == 0) {
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
  }
  
  // Final Result
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance,
                                 Rcpp::Named("Dev") = Dev,
                                 Rcpp::Named("Mon") = Mon,
                                 Rcpp::Named("thinned") = thinned));
}

// [[Rcpp::export]]
SEXP harm(Function Model, List Data, int Iterations, int Status,
          int Thinning, double Acceptance, NumericMatrix Dev,
          int LIV, NumericMatrix Mon, List Mo0, NumericMatrix thinned) {
  
  // Initial settings
  int t_iter = 0;
  double alpha = 0;
  Rcpp::List Mo1 = clone(Mo0);
  RNGScope scope;
  
  // Run MCMC algorithm
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 <<
        ",   Proposal: Multivariate,   LP: " <<
          floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
      t_iter = floor((iter) / Thinning) + 1;
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
    // Propose new values
    Rcpp::List          Mo0_ = clone(Mo0);
    Rcpp::NumericVector prop = HARproposal(Mo0_["parm"]);
    Rcpp::List          Mo1 = Model(prop, Data);
    // Accept/Reject
    double u = as<double>(runif(1));
    double LP0 = Mo0_["LP"];
    double LP1 = Mo1["LP"];
    alpha = exp(LP1 - LP0);
    if (u < alpha) {
      Mo0 = Mo1;
      Acceptance += 1.0 / LIV;
    }
    if ((iter + 1) % Thinning == 0) {
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
  }
  
  // Final Result
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance,
                                 Rcpp::Named("Dev") = Dev,
                                 Rcpp::Named("Mon") = Mon,
                                 Rcpp::Named("thinned") = thinned));
}

// [[Rcpp::export]]
SEXP sharm(Function Model, List Data, int Iterations, int Status,
           int Thinning, double Acceptance, NumericMatrix Dev, double h,
           int LIV, NumericMatrix Mon, List Mo0, NumericMatrix thinned) {
  
  // Initial settings
  int t_iter = 0;
  double alpha = 0;
  Rcpp::List Mo1 = clone(Mo0);
  RNGScope scope;
  NumericVector prop0 = as<Rcpp::NumericVector>(Mo0["parm"]);
  NumericVector gr0 = grad(Model, Data, as<Rcpp::NumericVector>(Mo0["parm"]), h);
  
  // Run MCMC algorithm
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 <<
        ",   Proposal: Multivariate,   LP: " <<
          floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
    }
    // Propose new values
    Rcpp::NumericVector prop = SHARproposal(as<Rcpp::NumericVector>(Mo0["parm"]), h, Data, Model, gr0);
    Rcpp::List Mo1 = Model(prop, Data);
    // Accept/Reject
    double u = as<double>(runif(1));
    double LP0 = Mo0["LP"];
    double LP1 = Mo1["LP"];
    alpha = exp(LP1 - LP0);
    if (u < alpha) {
        Mo0 = Mo1;
        gr0 = grad(Model, Data, prop, h);
        Acceptance += 1.0 / LIV;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
      t_iter = floor((iter) / Thinning) + 1;
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
    }
  }
  
  // Final Result
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance,
                                 Rcpp::Named("Dev") = Dev,
                                 Rcpp::Named("Mon") = Mon,
                                 Rcpp::Named("thinned") = thinned));
}

// ====================================THE END====================================