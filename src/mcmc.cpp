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

NumericMatrix Leapfrog(NumericVector par, NumericVector r, double epsilon,
                       Function Model, List Data, double h) {
  NumericVector r_hat = r + (epsilon / 2) * grad(Model, Data, par, h);
  NumericVector theta_hat = par + (epsilon * r_hat);
  NumericVector r_tilt = r_hat + (epsilon / 2) * grad(Model, Data, theta_hat, h);
  NumericMatrix Result(par.length(), 2);
  Result(_, 0) = theta_hat;
  Result(_, 1) = r_tilt;
  return Result;
}

arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  arma::mat m = Y * arma::chol(sigma);
  return m;
}

double FindReasonableEpsilon(NumericVector par, Function Model, List Data, double h) {
  // Initialize
  double        epsilon = 1.0;
  NumericVector mu(par.length());
  arma::mat     Sigma = Rcpp::as<arma::mat>(NumericMatrix::diag(par.length(), 1));
  NumericVector r = as<NumericVector>(wrap(mvrnormArma(1, mu, Sigma)));
  // Leapfrog step
  NumericMatrix LF = as<NumericMatrix>(Leapfrog(par, r, epsilon, Model, Data, h));
  NumericVector theta_tilt = as<NumericVector>(wrap(LF(_, 0)));
  NumericVector r_tilt = as<NumericVector>(wrap(LF(_, 1)));
  double        r_m = sum(r * r);
  double        rt_m = sum(r_tilt * r_tilt);
  double pU = exp(as<double>(as<List>(Model(theta_tilt, Data))["LP"]) - (.5 * rt_m));
  double pB = exp(as<double>(as<List>(Model(par, Data))["LP"]) - (.5 * r_m));
  // Establish the value of a
  double decision;
  if ((pU / pB) > .5) {
    decision = 1.0;
  }
  else {
    decision = 0.0;
  };
  double a = 2.0 * decision - 1.0;
  // Iterate
  while (pow((pU / pB), a) > pow(2, -a)) {
    epsilon = pow(2, a) * epsilon;
    // Leapfrog step
    LF = as<NumericMatrix>(Leapfrog(par, r, epsilon, Model, Data, h));
    theta_tilt = as<NumericVector>(wrap(LF(_, 0)));
    r_tilt = as<NumericVector>(wrap(LF(_, 1)));
    r_m = sum(r * r);
    rt_m = sum(r_tilt * r_tilt);
    pU = exp(as<double>(as<List>(Model(theta_tilt, Data))["LP"]) - (.5 * rt_m));
    pB = exp(as<double>(as<List>(Model(par, Data))["LP"]) - (.5 * r_m));
  }
  
  return epsilon;
}

List BuildTree(NumericVector theta, NumericVector r, double u, int v,
               int j, double epsilon, NumericVector theta_naught,
               double rn_m, Function Model, List Data, double h) {
  double delta_max = 1000.0;
  if (j == 0) {
    NumericMatrix LF = as<NumericMatrix>(Leapfrog(theta, r, ((double)v) * epsilon,
                                                  Model, Data, h));
    NumericVector theta_tilt = as<NumericVector>(wrap(LF(_, 0)));
    NumericVector r_tilt = as<NumericVector>(wrap(LF(_, 1)));
    double        rt_m = sum(r_tilt * r_tilt);
    int n_tilt;
    if (u <= exp(as<double>(as<List>(Model(theta_tilt, Data))["LP"]) - (0.5 * rt_m))) {
      n_tilt = 1;
    }
    else {
      n_tilt = 0;
    };
    int s_tilt;
    if (u < exp(delta_max + as<double>(as<List>(Model(theta_tilt, Data))["LP"]) - (0.5 * rt_m))) {
      s_tilt = 1;
    }
    else {
      s_tilt = 0;
    };
    double LPmin = as<double>(as<List>(Model(theta_tilt, Data))["LP"]) - (0.5 * rt_m) -
      as<double>(as<List>(Model(theta_naught, Data))["LP"]) - (0.5 * rn_m);
    double alpha_tilt = std::min(1.0, exp(LPmin));
    int    n_alpha = 1;
    
    return wrap(Rcpp::List::create(Rcpp::Named("theta_minus") = theta_tilt,
                                   Rcpp::Named("r_minus") = r_tilt,
                                   Rcpp::Named("theta_plus") = theta_tilt,
                                   Rcpp::Named("r_plus") = r_tilt,
                                   Rcpp::Named("theta_tilt") = theta_tilt,
                                   Rcpp::Named("n_tilt") = n_tilt,
                                   Rcpp::Named("s_tilt") = s_tilt,
                                   Rcpp::Named("alpha_tilt") = alpha_tilt,
                                   Rcpp::Named("n_alpha") = n_alpha));
  }
  else {
    List BT = BuildTree(theta, r, u, v, j - 1, epsilon, theta_naught, rn_m, Model, Data, h);
    NumericVector theta_minus = BT["theta_minus"];
    NumericVector r_minus = BT["r_minus"];
    NumericVector theta_plus = BT["theta_plus"];
    NumericVector r_plus = BT["r_plus"];
    NumericVector theta_tilt = BT["theta_tilt"];
    int n_tilt = BT["n_tilt"];
    int s_tilt = BT["s_tilt"];
    double alpha_tilt = BT["alpha_tilt"];
    int    n_alpha = BT["n_alpha"];
    if (s_tilt == 1) {
      NumericVector theta_Dtilt;
      int n_Dtilt;
      int s_Dtilt;
      double alpha_Dtilt;
      int    n_Dalpha;
      if (v == -1) {
        List BT = BuildTree(theta_minus, r_minus, u, v, j - 1, epsilon, theta_naught,
                            rn_m, Model, Data, h);
        NumericVector theta_minus = BT["theta_minus"];
        NumericVector r_minus = BT["r_minus"];
        NumericVector theta_Dtilt = BT["theta_tilt"];
        n_Dtilt = BT["n_tilt"];
        s_Dtilt = BT["s_tilt"];
        alpha_Dtilt = BT["alpha_tilt"];
        n_Dalpha = BT["n_alpha"];
      }
      else {
        List BT = BuildTree(theta_plus, r_plus, u, v, j - 1, epsilon, theta_naught,
                            rn_m, Model, Data, h);
        NumericVector theta_plus = BT["theta_plus"];
        NumericVector r_plus = BT["r_plus"];
        NumericVector theta_Dtilt = BT["theta_tilt"];
        n_Dtilt = BT["n_tilt"];
        s_Dtilt = BT["s_tilt"];
        alpha_Dtilt = BT["alpha_tilt"];
        n_Dalpha = BT["n_alpha"];
      }
      double prob = ((double)n_Dtilt) / (((double)n_tilt) + ((double)n_Dtilt));
      double coin = as<double>(runif(1));
      if (prob > coin) {
        theta_tilt = theta_Dtilt;
      }
      alpha_tilt += alpha_Dtilt;//alpha_tilt + alpha_Dtilt;
      n_alpha += n_Dalpha;//n_alpha + n_Dalpha;
      int rule1;
      if (sum((theta_plus - theta_minus) * r_minus) >= 0.0) {
        rule1 = 1;
      }
      else {
        rule1 = 0;
      }
      int rule2;
      if (sum((theta_plus - theta_minus) * r_plus) >= 0.0) {
        rule2 = 1;
      }
      else {
        rule2 = 0;
      }
      int rule;
      if ((rule1 + rule2) == 2) {
        rule = 1;
      }
      else {
        rule = 0;
      }
      s_tilt = s_Dtilt * rule;
      n_tilt += n_Dtilt;//n_tilt + n_Dtilt;
    }
    
    return wrap(Rcpp::List::create(Rcpp::Named("theta_minus") = theta_minus,
                                   Rcpp::Named("r_minus") = r_minus,
                                   Rcpp::Named("theta_plus") = theta_plus,
                                   Rcpp::Named("r_plus") = r_plus,
                                   Rcpp::Named("theta_tilt") = theta_tilt,
                                   Rcpp::Named("n_tilt") = n_tilt,
                                   Rcpp::Named("s_tilt") = s_tilt,
                                   Rcpp::Named("alpha_tilt") = alpha_tilt,
                                   Rcpp::Named("n_alpha") = n_alpha));
  }
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

Rcpp::NumericVector SHARproposal(NumericVector par, double h, List Data, Function Model, double max_d) {
  
  Rcpp::NumericVector d = runif(par.length(), 0, max_d);
  Rcpp::NumericVector gr = grad(Model, Data, par, h);
  Rcpp::NumericVector newgr = gr / sqrt(sum(gr * gr));
  Rcpp::NumericVector prop = HARproposal(par + (d * newgr));
  
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
  double max_d = FindReasonableEpsilon(Mo0["parm"], Model, Data, h);

  
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
    Rcpp::NumericVector prop = SHARproposal(Mo0_["parm"], h, Data, Model, max_d);
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
SEXP nutsda(Function Model, List Data, int Iterations, int Status, double h,
            int Thinning, double Acceptance, NumericMatrix Dev, int M_adap,
            int LIV, NumericMatrix Mon, List Mo0, NumericMatrix thinned) {
  
  // Initial settings
  int    t_iter = 0;
  double epsilon_0 = FindReasonableEpsilon(Mo0["parm"], Model, Data, h);
  NumericVector epsilon_m(Iterations);
  epsilon_m[0] = epsilon_0;
  double mu = log(10 * epsilon_0);
  NumericVector epsilon_bar(Iterations);
  epsilon_bar[0] = 1.0;
  NumericVector H_bar(Iterations);
  H_bar[0] = 0.0;
  double gamma = 0.05;
  double t_0 = 10.0;
  double kappa = .75;
  NumericVector par_0 = as<NumericVector>(Mo0["parm"]);
  NumericVector expectation(par_0.length());
  arma::mat     Sigma = Rcpp::as<arma::mat>(NumericMatrix::diag(par_0.length(), 1));
  Rcpp::List Mo1 = clone(Mo0);
  RNGScope scope;
  
  // Run MCMC algorithm
  for (int iter = 0; iter < Iterations; iter++) {
    // Print Status
    if ((iter + 1) % Status == 0) {
      Rcpp::Rcout << "Iteration: " << iter + 1 <<
        ",   Proposal: Multivariate,   LP: " <<
          floor(as<double>(Mo1["LP"]) * 100) / 100 << std::endl;
    }
    // Save Thinned Samples
    if ((iter + 1) % Thinning == 0) {
      t_iter = floor((iter) / Thinning) + 1;
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo1["parm"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo1["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo1["Monitor"]);
    }
    // Propose new values
    //NumericVector r_0 = as<NumericVector>(rnorm(par_0.length(), 0, 1));
    NumericVector r_0 = as<NumericVector>(wrap(mvrnormArma(1, expectation, Sigma)));
    // Build Tree for proposal
    double        rn_m = sum(r_0 * r_0);
    double        m_1 = as<double>(Mo1["LP"]);
    double        u = as<double>(runif(1)) * exp(m_1 - (.5 * rn_m));
    NumericVector theta_tilt;
    NumericVector theta_minus = Mo1["parm"];
    NumericVector theta_plus = Mo1["parm"];
    NumericVector r_minus = r_0;
    NumericVector r_plus = r_0;
    int j = 0;
    NumericVector theta_m = Mo1["parm"];
    int n = 1;
    int s = 1;
    int n_tilt = 1;
    int s_tilt = 1;
    double alpha;
    int    n_alpha;
    while (s == 1) {
      double coin = as<double>(runif(1));
      int v;
      if (coin >= .5) {
        v = 1;
      }
      else {
        v = -1;
      }
      List BT;
      if (v == -1) {
        BT = BuildTree(theta_minus, r_minus, u, v, j, epsilon_m[iter],
                       Mo1["parm"], rn_m, Model, Data, h);
        theta_minus = BT["theta_minus"];
        r_minus = BT["r_minus"];
        theta_tilt = BT["theta_tilt"];
        n_tilt = BT["n_tilt"];
        s_tilt = BT["s_tilt"];
        alpha = BT["alpha_tilt"];
        n_alpha = BT["n_alpha"];
      }
      else {
        BT = BuildTree(theta_plus, r_plus, u, v, j, epsilon_m[iter],
                       Mo1["parm"], rn_m, Model, Data, h);
        theta_plus = BT["theta_plus"];
        r_plus = BT["r_plus"];
        theta_tilt = BT["theta_tilt"];
        n_tilt = BT["n_tilt"];
        s_tilt = BT["s_tilt"];
        alpha = BT["alpha_tilt"];
        n_alpha = BT["n_alpha"];
      }
      if (s_tilt == 1) {
        double prob = std::min<double>(1.0, ((double)n_tilt / (double)n));
        double dec = as<double>(runif(1));
        if (prob > dec) {
          theta_m = theta_tilt;
        }
      }
      n += n_tilt;
      int rule1;
      if (sum((theta_plus - theta_minus) * r_minus) >= 0.0) {
        rule1 = 1;
      }
      else {
        rule1 = 0;
      }
      int rule2;
      if (sum((theta_plus - theta_minus) * r_plus) >= 0.0) {
        rule2 = 1;
      }
      else {
        rule2 = 0;
      }
      int rule;
      if ((rule1 + rule2) == 2) {
        rule = 1;
      }
      else {
        rule = 0;
      }
      s = s_tilt * rule;
      j += 1;
    }
    // Accept/Reject and Update
    Mo1 = Model(theta_m, Data);
    if (iter > 0) {
      if (iter <= M_adap) {
        H_bar[iter] = ((1.0 - (1.0 / (((double)iter) + t_0))) * H_bar[iter - 1]) +
          ((1.0 / (((double)iter) + t_0)) * (Acceptance - (alpha / ((double)n_alpha))));
        epsilon_m[iter] = exp(mu - (std::sqrt((double)iter) / gamma) * H_bar[iter]);
        epsilon_bar[iter] = exp(pow((double)iter, -kappa) * log(epsilon_m[iter]) +
          (1 - pow((double)iter, -kappa)) * log(epsilon_bar[iter - 1]));
      }
      else {
        epsilon_m[iter] = epsilon_m[iter - 1];
        epsilon_bar[iter] = epsilon_m[iter];
      }
    }
    if ((iter + 1) % Thinning == 0) {
      thinned(t_iter, _) = as<Rcpp::NumericVector>(Mo1["parm"]);
      Dev(t_iter, _) = as<Rcpp::NumericVector>(Mo1["Dev"]);
      Mon(t_iter, _) = as<Rcpp::NumericVector>(Mo1["Monitor"]);
    }
  }
  
  // Final Result
  return wrap(Rcpp::List::create(Rcpp::Named("Acceptance") = Acceptance,
                                 Rcpp::Named("Dev") = Dev,
                                 Rcpp::Named("Mon") = Mon,
                                 Rcpp::Named("thinned") = thinned));
}

// ====================================THE END====================================