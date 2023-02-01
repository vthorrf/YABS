// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(unwindProtect)]]

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

Rcpp::NumericMatrix covMat(const Rcpp::NumericMatrix &X) {
  arma::mat X_arma = Rcpp::as<arma::mat>(X);
  int n = X_arma.n_rows;
  arma::mat centered = X_arma - arma::repmat(arma::mean(X_arma, 0), n, 1);
  arma::mat covariance = (centered.t() * centered) / double(n - 1);
  return Rcpp::wrap(covariance);
}

Rcpp::NumericVector RWproposal(NumericVector par, NumericMatrix epsilon) {
  
  NumericVector mu(par.length());
  arma::mat     Sigma = Rcpp::as<arma::mat>( epsilon );
  Rcpp::NumericVector d = as<NumericVector>(wrap(mvrnormArma(1, mu, Sigma)));
  Rcpp::NumericVector prop = par + d;
  
  return prop;
}

Rcpp::NumericVector Bproposal(NumericVector par, NumericVector gr, NumericMatrix epsilon) {
  
  NumericVector mu(par.length());
  arma::mat     Sigma = Rcpp::as<arma::mat>( epsilon );
  Rcpp::NumericVector d = as<NumericVector>(wrap(mvrnormArma(1, mu, Sigma)));
  Rcpp::NumericVector Pr = 1 / (1 + exp(-(d * gr)));
  Rcpp::NumericVector u = runif(par.length());
  Rcpp::NumericVector B(par.length());
  for (int iter = 0; iter < par.length(); iter++) {
      if (u[iter] < Pr[iter]) {
          B[iter] = 1;
      }
      else {
          B[iter] = -1;
      }
  }
  Rcpp::NumericVector prop = par + (B * d);

  return prop;
}

// ====================================MCMC Algorithms====================================

// [[Rcpp::export]]
SEXP harmwg(Function Model, List Data, int Iterations, int Status,
            int Thinning, double ACC, NumericMatrix DevianceMat,
            int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples,
            NumericMatrix PPD, int Adapt, NumericMatrix Sigma) {
  
  // Initial settings
  double Acceptance = ACC;
  Rcpp::NumericMatrix Dev = clone(DevianceMat);
  Rcpp::NumericMatrix Mon = clone(Monitor);
  Rcpp::List Mo1 = clone(Mo0);
  Rcpp::NumericMatrix thinned = clone(samples);
  Rcpp::NumericMatrix postpred = clone(PPD);
  Rcpp::NumericMatrix epsilon = clone(Sigma);
  RNGScope scope;
  int t_iter = 0, fins = 0, mcols = Mon.ncol();
  double alpha = .0;
  
  // Run MCMC algorithm
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 <<
        ",   Proposal: Componentwise,   LP: " <<
          floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
    }
    // Random-Scan Componentwise Estimation
    Rcpp::NumericVector u = runif(LIV),      // Acceptance threshold for each parameter
                        z = rnorm(LIV);      // New Proposed value
    Rcpp::IntegerVector LIVseq = Rcpp::Range(0, LIV - 1), // Indexes to sample
                        s = Rcpp::RcppArmadillo::sample(LIVseq, LIV, false, NumericVector::create());  // Sample of indexes
    // Propose and evaluate new values per parameter
    for (int j = 0; j < LIV; j++) {
      // Propose new parameter values
      Rcpp::List Mo0_ = clone(Mo0);
      Rcpp::NumericVector prop = Mo0_["parm"];
      Rcpp::NumericVector prop1 = RWproposal(prop, epsilon);
      prop[s[j]] = prop1[s[j]];
      // Log-Posterior of the proposed state
      Rcpp::List Mo1 = Model(prop, Data);
      fins = ::R_finite(Mo1["LP"]) + ::R_finite(Mo1["Dev"]);
      for (int m = 0; m < mcols; m++) {
        fins += ::R_finite(as<Rcpp::NumericVector>(Mo1["LP"])[m]);
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
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
        t_iter = floor((iter + 1) / Thinning) - 1;
        thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
        postpred(t_iter, _) = as<Rcpp::NumericVector>(Mo0["yhat"]);
        Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
        Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
        if((iter < Adapt) & ((t_iter + 1) > 20)) {
          epsilon = covMat( thinned(Range(0,t_iter), _) );
        }
    }
  }
  
  // Final Result
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance / Iterations,
                                 Rcpp::Named("Dev") = Dev,
                                 Rcpp::Named("Mon") = Mon,
                                 Rcpp::Named("thinned") = thinned,
                                 Rcpp::Named("postpred") = postpred));
}

// [[Rcpp::export]]
SEXP harm(Function Model, List Data, int Iterations, int Status,
          int Thinning, double ACC, NumericMatrix DevianceMat,
          int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples,
          NumericMatrix PPD, int Adapt, NumericMatrix Sigma) {
  
  // Initial settings
  int t_iter = 0;
  double alpha = .0;
  double Acceptance = ACC;
  Rcpp::NumericMatrix Dev = clone(DevianceMat);
  Rcpp::NumericMatrix Mon = clone(Monitor);
  Rcpp::List Mo1 = clone(Mo0);
  Rcpp::NumericMatrix thinned = clone(samples);
  Rcpp::NumericMatrix postpred = clone(PPD);
  Rcpp::NumericMatrix epsilon = clone(Sigma);
  RNGScope scope;
  
  // Run MCMC algorithm
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 <<
        ",   Proposal: Multivariate,   LP: " <<
          floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
    }
    // Propose new values
    Rcpp::List          Mo0_ = clone(Mo0);
    Rcpp::NumericVector prop = RWproposal(Mo0_["parm"], epsilon);
    Rcpp::List          Mo1 = Model(prop, Data);
    // Accept/Reject
    double u = as<double>(runif(1));
    double LP0 = Mo0_["LP"];
    double LP1 = Mo1["LP"];
    alpha = exp(LP1 - LP0);
    if (u < alpha) {
      Mo0 = Mo1;
      Acceptance += 1.0;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
        t_iter = floor((iter + 1) / Thinning) - 1;
        thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
        postpred(t_iter, _) = as<Rcpp::NumericVector>(Mo0["yhat"]);
        Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
        Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
        if((iter < Adapt) & ((t_iter + 1) > 20)) {
          epsilon = covMat( thinned(Range(0,t_iter), _) );
        }
    }
  }
  
  // Final Result
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance / Iterations,
                                 Rcpp::Named("Dev") = Dev,
                                 Rcpp::Named("Mon") = Mon,
                                 Rcpp::Named("thinned") = thinned,
                                 Rcpp::Named("postpred") = postpred));
}

// [[Rcpp::export]]
SEXP gcharm(Function Model, List Data, int Iterations, int Status,
            int Thinning, double ACC, NumericMatrix DevianceMat, double h,
            int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples,
            NumericMatrix PPD, int Adapt, NumericMatrix Sigma) {
  
  // Initial settings
  int t_iter = 0;
  double alpha = 0;
  double Acceptance = ACC;
  Rcpp::NumericMatrix Dev = clone(DevianceMat);
  Rcpp::NumericMatrix Mon = clone(Monitor);
  Rcpp::List Mo1 = clone(Mo0);
  Rcpp::NumericMatrix thinned = clone(samples);
  Rcpp::NumericMatrix postpred = clone(PPD);
  Rcpp::NumericMatrix epsilon = clone(Sigma);
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
    Rcpp::NumericVector prop = Bproposal(as<Rcpp::NumericVector>(Mo0["parm"]), gr0, epsilon);
    Rcpp::List Mo1 = Model(prop, Data);
    // Accept/Reject
    double u = as<double>(runif(1));
    double LP0 = Mo0["LP"];
    double LP1 = Mo1["LP"];
    alpha = 1 / (exp(LP0 - LP1) + 1);
    if (u < alpha) {
        Mo0 = Mo1;
        gr0 = grad(Model, Data, prop, h);
        Acceptance += 1.0;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
      t_iter = floor((iter + 1) / Thinning) - 1;
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
      postpred(t_iter, _) = as<Rcpp::NumericVector>(Mo0["yhat"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
      if((iter < Adapt) & ((t_iter + 1) > 20)) {
        epsilon = covMat( thinned(Range(0,t_iter), _) );
      }
    }
  }
  
  // Final Result
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance / Iterations,
                                 Rcpp::Named("Dev") = Dev,
                                 Rcpp::Named("Mon") = Mon,
                                 Rcpp::Named("thinned") = thinned,
                                 Rcpp::Named("postpred") = postpred));
}

// [[Rcpp::export]]
SEXP ohss(Function Model, List Data, int Iterations, int Status,
          int Thinning, double ACC, NumericMatrix DevianceMat,
          int LIV, NumericMatrix Monitor, List Mo0, NumericMatrix samples,
          NumericMatrix PPD, int Adapt) {

    // Initial settings
    double Acceptance = ACC;
    Rcpp::NumericMatrix Dev = clone(DevianceMat);
    Rcpp::NumericMatrix Mon = clone(Monitor);
    Rcpp::List Mo1 = clone(Mo0);
    Rcpp::NumericMatrix thinned = clone(samples);
    Rcpp::NumericMatrix postpred = clone(PPD);
    Rcpp::NumericMatrix post = thinned;
    post(0, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
    int t_iter = 0;
    double y_slice = Mo0["LP"];
    Rcpp::NumericVector L(LIV);
    Rcpp::NumericVector U(LIV);
    Rcpp::NumericVector wt(LIV);
    Rcpp::NumericVector v(LIV);
    double w = 0.05;
    Rcpp::NumericMatrix VarCov = NumericMatrix::diag(LIV, 1);
    Rcpp::NumericMatrix VarCov2 = VarCov;
    int decomp_freq = max(IntegerVector::create((int)floor(Iterations/Thinning/100), 10));
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, as<arma::mat>(wrap(VarCov)));
    Rcpp::NumericMatrix S_eig = as<Rcpp::NumericMatrix>(wrap(eigvec));
    Rcpp::NumericVector V_eig = as<Rcpp::NumericVector>(wrap(eigval)).sort(true);
    Rcpp::NumericMatrix DiagCovar((int)floor(Iterations/Thinning) + 1, LIV);
    for (int i = 0; i < LIV; i++) {
        DiagCovar(0,i) = S_eig(i, i);
    }
    double tuning = 1.0;
    double edge_scale = 5.0;
    RNGScope scope;

    // Run MCMC algorithm
    for (int iter = 0; iter < Iterations; iter++) {
        // Print Status
        if ((iter + 1) % Status == 0) {
            Rcpp::Rcout << "Iteration: " << iter + 1 <<
                ",   Proposal: Multivariate,   LP: " <<
                floor(as<double>(Mo0["LP"]) * 100) / 100 << std::endl;
        }
        // Eigenvectors of the Sample Covariance Matrix
        if ( ((iter % decomp_freq) == 0) && (iter > 1) && (iter < Adapt) ) {
            arma::mat temp = as<arma::mat>(wrap(post(Range(0, iter - 1), _)));
            VarCov2 = Rcpp::NumericMatrix(wrap(temp.t() * temp)) / (iter - 1);
            arma::vec eigval;
            arma::mat eigvec;
            arma::eig_sym(eigval, eigvec, as<arma::mat>(wrap(VarCov2)));
        }
        // Hypercube or Eigenvector
        double ru = as<double>(runif(1));
        if (ru < w) {
            V_eig = rep(tuning, LIV);
            S_eig = NumericMatrix::diag(LIV, 1);
        }
        else {
            V_eig = as<Rcpp::NumericVector>(wrap(eigval)).sort(true);
            S_eig = as<Rcpp::NumericMatrix>(wrap(eigvec));
        }
        // Slice Interval
        y_slice = as<double>(Mo0["LP"]) - as<double>(rexp(1, 1.0));
        L = -1.0 * runif(LIV);
        U = L + 1.0;
        // Rejection Sampling
        int condition = 0;
        while (condition != 1) {
            for(int jter = 0; jter < LIV; jter++) {
              wt[jter] = as<double>(runif(1, L[jter], U[jter]));
            }
            v = as<arma::mat>(wrap(S_eig)) * (edge_scale * as<arma::vec>(wrap(wt)) % as<arma::vec>(wrap(V_eig)));
            NumericVector prop = as<Rcpp::NumericVector>(Mo0["parm"]) + as<Rcpp::NumericVector>(wrap(v));
            Mo1 = Model(prop, Data);
            int ALL = 0;
            for (int jter = 0; jter < LIV; jter++) {
                if (abs(wt[jter]) < 1e-100) {
                    ALL += 1;
                }
            }
            if (as<double>(wrap(Mo1["LP"])) >= y_slice) {
                condition = 1;
                break;
            }
            else if (ALL == LIV) {
                condition = 1;
                Mo1 = Mo0;
                break;
            }
            for (int jter = 0; jter < LIV; jter++) {
                if (wt[jter] < 0.0) {
                    L[jter] = wt[jter];
                }
                if (wt[jter] > 0.0) {
                    U[jter] = wt[jter];
                }
            }
        }
        Mo0 = Mo1;
        if (iter < Adapt) {
            post(iter,_) = as<Rcpp::NumericVector>(Mo0["parm"]);
            VarCov = VarCov2;
        }
        // Save Thinned Samples
        Acceptance += 1.0;
        if ((iter + 1) % Thinning == 0) {
            t_iter = floor((iter + 1) / Thinning) - 1;
            thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo0["parm"]);
            postpred(t_iter, _) = as<Rcpp::NumericVector>(Mo0["yhat"]);
            Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Dev"]);
            Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo0["Monitor"]);
        }
    }

    // Final Result
    return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance / Iterations,
                                   Rcpp::Named("Dev") = Dev,
                                   Rcpp::Named("Mon") = Mon,
                                   Rcpp::Named("thinned") = thinned,
                                   Rcpp::Named("postpred") = postpred));
}

// ====================================THE END====================================