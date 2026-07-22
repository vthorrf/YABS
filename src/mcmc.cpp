// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using namespace arma;

// ====================================Supporting functions====================================
arma::vec grad(Function Model, List Data, arma::vec par, double h = 1e-6) {

    int m = par.n_elem;
    arma::vec out(m);

    double f_x = as<double>(as<List>(Model(wrap(par), Data))["LP"]);

    for (int i = 0; i < m; i++) {
        arma::vec temp = par;
        temp[i] += h;

        double f1 = as<double>(as<List>(Model(wrap(temp), Data))["LP"]);
        out[i] = (f1 - f_x) / h;

        if (!std::isfinite(out[i])) out[i] = 0.0;
    }

    return out;
}

arma::vec gradC(Function Model, List Data, const arma::vec& par, double h = 1e-6) {

    int m = par.n_elem;
    arma::vec out(m, fill::zeros);

    for (int i = 0; i < m; i++) {

        arma::vec forward = par;
        arma::vec backward = par;

        forward[i] += h;
        backward[i] -= h;

        double f_plus = as<double>(as<List>(Model(wrap(forward), Data))["LP"]);
        double f_minus = as<double>(as<List>(Model(wrap(backward), Data))["LP"]);

        if (!std::isfinite(f_plus))  f_plus = -1e300;
        if (!std::isfinite(f_minus)) f_minus = -1e300;

        out[i] = (f_plus - f_minus) / (2.0 * h);

        if (!std::isfinite(out[i])) out[i] = 0.0;
    }

    return out;
}

// [[Rcpp::export]]
arma::vec gradN(Function Model, List Data, arma::vec par,
                double h = 1e-6, int order = 1) {

    int m = par.n_elem;
    arma::vec out(m, fill::zeros);

    // Precompute base evaluation
    double f0 = as<double>(as<List>(Model(wrap(par), Data))["LP"]);

    // Precompute binomial coefficients
    arma::vec binom(order + 1, fill::zeros);
    binom[0] = 1.0;
    for (int i = 1; i <= order; i++) {
        binom[i] = binom[i - 1] * (order - i + 1) / i;
    }

    for (int i = 0; i < m; i++) {

        double acc = 0.0;

        for (int j = 0; j <= order; j++) {

            arma::vec temp = par;
            temp[i] += j * h;

            double fj;

            if (j == 0) {
                fj = f0;
            }
            else {
                fj = as<double>(as<List>(Model(wrap(temp), Data))["LP"]);
            }

            if (!std::isfinite(fj)) fj = 0.0;

            // finite difference weight
            double sign = ((order - j) % 2 == 0) ? 1.0 : -1.0;
            acc += sign * binom[j] * fj;
        }

        out[i] = acc / std::pow(h, order);

        if (!std::isfinite(out[i])) out[i] = 0.0;
    }

    return out;
}

arma::mat mvrnormArma(int n, const arma::mat& Sigma) {
    int p = Sigma.n_cols;
    arma::mat Y = arma::randn(n, p);
    arma::mat C = arma::chol(Sigma);
    return Y * C;
}

arma::vec RWproposal(const arma::vec& par, const arma::mat& Sigma) {
    arma::vec d = mvrnormArma(1, Sigma).row(0).t();
    return par + d;
}

arma::vec Bproposal(const arma::vec& par, const arma::vec& gr,
                    const arma::mat& Sigma) {

    arma::vec d = mvrnormArma(1, Sigma).row(0).t();
    arma::vec Pr = 1.0 / (1.0 + arma::exp(-(d % gr)));

    arma::vec u = arma::randu<arma::vec>(par.n_elem);
    arma::vec B = arma::ones<arma::vec>(par.n_elem);
    B.elem(find(u >= Pr)).fill(-1.0);

    return par + (B % d);
}

arma::mat sigmaEstimate(const arma::mat& X, const arma::mat& epsilon,
                        int ind, const arma::vec& sigma, bool accept) {

    double pi = datum::pi;
    double p_star = 0.40;
    double alpha = 0.8416212;

    double c1 = (1.0 - (1.0 / sigma.n_elem));
    double c2 = (std::sqrt(2.0 * pi) * std::exp(alpha * alpha / 2.0)) / (2.0 * alpha);
    double c3 = 1.0 / (sigma.n_elem * p_star * (1 - p_star));
    
    double factor = c1 * c2 + c3;
    
    arma::vec stepLength = sigma * factor;

    arma::vec sigma_new = sigma;

    if (accept) {
        sigma_new += (stepLength * (1 - p_star)) /
            std::max(200.0, double(ind / sigma.n_elem));
    }
    else {
        sigma_new -= (stepLength * p_star) /
            std::max(200.0, double(ind / sigma.n_elem));
    }

    arma::mat Id = diagmat(square(sigma_new) / ind);

    if (ind <= 100) {
        return chol(epsilon * epsilon.t() + Id).t();
    }
    else {
        arma::mat centered = X.each_row() - mean(X, 0);
        arma::mat cov = (centered.t() * centered) / (X.n_rows - 1);
        return chol(cov + Id).t();
    }
}

arma::mat Leapfrog(const arma::vec& theta, const arma::vec& r,
                   double epsilon, Function Model, List Data,
                   double h, const arma::vec& grad_theta) {

    arma::vec r_half = r + 0.5 * epsilon * grad_theta;
    arma::vec theta_new = theta + epsilon * r_half;

    arma::vec g2 = gradC(Model, Data, theta_new, h);
    arma::vec r_new = r_half + 0.5 * epsilon * g2;

    arma::mat out(theta.n_elem, 3);
    out.col(0) = theta_new;
    out.col(1) = r_new;
    out.col(2) = g2;

    return out;
}

List BuildTree(const arma::vec& theta, const arma::vec& r,
               const arma::vec& grad0, double logu, int v, int j,
               double epsilon, const arma::vec& theta0,
               double joint0, Function Model, List Data, double h) {

    const double delta_max = 1000.0;

    // ================= BASE CASE =================
    if (j == 0) {

        arma::mat LF = Leapfrog(theta, r, v * epsilon, Model, Data, h, grad0);

        arma::vec theta1 = LF.col(0);
        arma::vec r1 = LF.col(1);
        arma::vec g1 = LF.col(2);

        List Mo1 = Model(theta1, Data);
        double lp1 = as<double>(Mo1["LP"]);
        double joint = lp1 - 0.5 * dot(r1, r1);

        int n = (logu <= joint) ? 1 : 0;
        int s = (logu < joint + delta_max) ? 1 : 0;

        double alpha = std::exp(std::min(0.0, joint - joint0));
        int n_alpha = 1;

        return List::create(
            Named("theta_minus") = theta1,
            Named("r_minus") = r1,
            Named("grad_minus") = g1,
            Named("theta_plus") = theta1,
            Named("r_plus") = r1,
            Named("grad_plus") = g1,
            Named("theta") = theta1,
            Named("n") = n,
            Named("s") = s,
            Named("alpha") = alpha,
            Named("n_alpha") = n_alpha
        );
    }

    // ================= RECURSION =================
    List left = BuildTree(theta, r, grad0, logu, v, j - 1, epsilon,
        theta0, joint0, Model, Data, h);

    arma::vec theta_minus = as<arma::vec>(left["theta_minus"]);
    arma::vec r_minus = as<arma::vec>(left["r_minus"]);
    arma::vec grad_minus = as<arma::vec>(left["grad_minus"]);
    arma::vec theta_plus = as<arma::vec>(left["theta_plus"]);
    arma::vec r_plus = as<arma::vec>(left["r_plus"]);
    arma::vec grad_plus = as<arma::vec>(left["grad_plus"]);
    arma::vec theta_prop = as<arma::vec>(left["theta"]);

    int    n = as<int>(left["n"]);
    int    s = as<int>(left["s"]);
    double alpha = as<double>(left["alpha"]);
    int    n_alpha = as<int>(left["n_alpha"]);

    if (s == 1) {

        List right;

        if (v == -1) {
            right = BuildTree(theta_minus, r_minus, grad_minus,
                logu, v, j - 1, epsilon,
                theta0, joint0, Model, Data, h);

            theta_minus = as<arma::vec>(right["theta_minus"]);
            r_minus = as<arma::vec>(right["r_minus"]);
            grad_minus = as<arma::vec>(right["grad_minus"]);
        }
        else {
            right = BuildTree(theta_plus, r_plus, grad_plus,
                logu, v, j - 1, epsilon,
                theta0, joint0, Model, Data, h);

            theta_plus = as<arma::vec>(right["theta_plus"]);
            r_plus = as<arma::vec>(right["r_plus"]);
            grad_plus = as<arma::vec>(right["grad_plus"]);
        }

        arma::vec theta2 = as<arma::vec>(right["theta"]);
        int    n2 = as<int>(right["n"]);
        int    s2 = as<int>(right["s"]);
        double alpha2 = as<double>(right["alpha"]);
        int    n_alpha2 = as<int>(right["n_alpha"]);

        if ((n + n2) > 0 && R::runif(0, 1) < (double)n2 / (n + n2)) {
            theta_prop = theta2;
        }

        n += n2;

        double dot1 = dot(theta_plus - theta_minus, r_minus);
        double dot2 = dot(theta_plus - theta_minus, r_plus);
        s = s2 && (dot1 >= 0.0) && (dot2 >= 0.0);

        alpha += alpha2;
        n_alpha += n_alpha2;
    }

    return List::create(
        Named("theta_minus") = theta_minus,
        Named("r_minus") = r_minus,
        Named("grad_minus") = grad_minus,
        Named("theta_plus") = theta_plus,
        Named("r_plus") = r_plus,
        Named("grad_plus") = grad_plus,
        Named("theta") = theta_prop,
        Named("n") = n,
        Named("s") = s,
        Named("alpha") = alpha,
        Named("n_alpha") = n_alpha
    );
}

double FindReasonableEpsilon(const arma::vec& theta, Function Model,
                             List Data, double h = 1e-6) {

    int dim = theta.n_elem;

    double epsilon = 1.0;
    arma::vec r = arma::randn(dim);

    List Mo0 = Model(theta, Data);
    double joint0 = as<double>(Mo0["LP"]) - 0.5 * dot(r, r);

    arma::vec g0 = gradC(Model, Data, theta, h);
    arma::mat LF = Leapfrog(theta, r, epsilon, Model, Data, h, g0);

    arma::vec theta_new = LF.col(0);
    arma::vec r_new = LF.col(1);

    List Mo1 = Model(theta_new, Data);
    double joint1 = as<double>(Mo1["LP"]) - 0.5 * dot(r_new, r_new);

    double log_accept = joint1 - joint0;
    int a = (log_accept > log(0.5)) ? 1 : -1;

    while (true) {

        epsilon = (a == 1) ? epsilon * 2.0 : epsilon * 0.5;

        LF = Leapfrog(theta, r, epsilon, Model, Data, h, g0);

        theta_new = LF.col(0);
        r_new = LF.col(1);

        Mo1 = Model(theta_new, Data);
        joint1 = as<double>(Mo1["LP"]) - 0.5 * dot(r_new, r_new);

        log_accept = joint1 - joint0;

        if (a == 1 && log_accept <= log(0.5)) break;
        if (a == -1 && log_accept >= log(0.5)) break;

        if (epsilon < 1e-10 || epsilon > 1e10) break;
    }

    return epsilon;
}

// ====================================MCMC Algorithms====================================

// [[Rcpp::export]]
List harmwg(Function Model, List Data, int Iterations, int Status,
            arma::vec InitialValues, int Thinning, double ACC,
            arma::mat DevianceMat, int LIV, arma::mat Monitor,
            List Mo0, arma::mat samples, arma::mat PPD,
            int Adapt, arma::mat Sigma) {

    RNGScope scope;

    arma::vec auxiliary = InitialValues;
    arma::mat epsilon = Sigma;

    arma::mat thinned = samples;
    arma::mat postpred = PPD;
    arma::mat Dev = DevianceMat;
    arma::mat Mon = Monitor;

    double Acceptance = ACC;
    int t_iter = 0;

    for (int iter = 0; iter < Iterations; iter++) {

        if ((iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                << ", LP: " << as<double>(Mo0["LP"]) << std::endl;
        }

        arma::vec prop = auxiliary;

        // random scan order
        arma::uvec idx = randperm(LIV);

        for (int j = 0; j < LIV; j++) {

            arma::vec prop_try = prop;

            arma::vec step = mvrnormArma(1, epsilon).row(0).t();
            prop_try[idx[j]] += step[idx[j]];

            List Mo1 = Model(wrap(prop_try), Data);

            double LP0 = as<double>(Mo0["LP"]);
            double LP1 = as<double>(Mo1["LP"]);

            if (!std::isfinite(LP1)) LP1 = -1e300;

            double alpha = std::exp(LP1 - LP0);
            double u = R::runif(0, 1);

            if (u < alpha) {
                prop = prop_try;
                Mo0 = Mo1;
                Acceptance += 1.0 / LIV;
            }
        }

        auxiliary = prop;

        if ((iter + 1) % Thinning == 0) {
            t_iter = (iter + 1) / Thinning - 1;

            thinned.row(t_iter) = as<arma::vec>(Mo0["parm"]).t();
            postpred.row(t_iter) = as<arma::vec>(Mo0["yhat"]).t();
            Dev.row(t_iter) = as<arma::vec>(Mo0["Dev"]).t();
            Mon.row(t_iter) = as<arma::vec>(Mo0["Monitor"]).t();

            if (iter < Adapt) {
                arma::vec sigma = sqrt(diagvec(epsilon * epsilon.t()));
                bool accept = (Acceptance > 0.5);

                epsilon = sigmaEstimate(
                    thinned.rows(0, t_iter),
                    epsilon,
                    t_iter + 1,
                    sigma,
                    accept
                );
            }
        }
    }

    return List::create(
        Named("Acceptance") = Acceptance / Iterations,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("thinned") = thinned,
        Named("postpred") = postpred
    );
}

// [[Rcpp::export]]
List harm(Function Model, List Data, int Iterations, int Status,
          arma::vec InitialValues, int Thinning, double ACC,
          arma::mat DevianceMat, int LIV, arma::mat Monitor,
          List Mo0, arma::mat samples, arma::mat PPD, int Adapt,
          arma::mat Sigma) {

    RNGScope scope;

    arma::vec auxiliary = InitialValues;
    arma::mat epsilon = Sigma;

    arma::mat thinned = samples;
    arma::mat postpred = PPD;
    arma::mat Dev = DevianceMat;
    arma::mat Mon = Monitor;

    double Acceptance = ACC;
    int t_iter = 0;

    for (int iter = 0; iter < Iterations; iter++) {

        if ((iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                << ", LP: " << as<double>(Mo0["LP"]) << std::endl;
        }

        arma::vec prop = RWproposal(auxiliary, epsilon);

        List Mo1 = Model(wrap(prop), Data);

        double LP0 = as<double>(Mo0["LP"]);
        double LP1 = as<double>(Mo1["LP"]);

        double alpha = std::exp(LP1 - LP0);
        double u = R::runif(0, 1);

        if (u < alpha) {
            Mo0 = Mo1;
            auxiliary = prop;
            Acceptance += 1.0;
        }

        if ((iter + 1) % Thinning == 0) {
            t_iter = (iter + 1) / Thinning - 1;

            thinned.row(t_iter) = as<arma::vec>(Mo0["parm"]).t();
            postpred.row(t_iter) = as<arma::vec>(Mo0["yhat"]).t();
            Dev.row(t_iter) = as<arma::vec>(Mo0["Dev"]).t();
            Mon.row(t_iter) = as<arma::vec>(Mo0["Monitor"]).t();

            if (iter < Adapt) {
                arma::vec sigma = sqrt(diagvec(epsilon * epsilon.t()));
                bool accept = (u < alpha);

                epsilon = sigmaEstimate(
                    thinned.rows(0, t_iter),
                    epsilon,
                    t_iter + 1,
                    sigma,
                    accept
                );
            }
        }
    }

    return List::create(
        Named("Acceptance") = Acceptance / Iterations,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("thinned") = thinned,
        Named("postpred") = postpred
    );
}

// [[Rcpp::export]]
List gcharm(Function Model, List Data, int Iterations, int Status,
            arma::vec InitialValues, int Thinning, double ACC,
            arma::mat DevianceMat, double h, int LIV,
            arma::mat Monitor, List Mo0, arma::mat samples,
            arma::mat PPD, int Adapt, arma::mat Sigma) {

    RNGScope scope;

    arma::vec auxiliary = InitialValues;
    arma::mat epsilon = Sigma;

    arma::vec gr0 = grad(Model, Data, auxiliary, h);

    arma::mat thinned = samples;
    arma::mat postpred = PPD;
    arma::mat Dev = DevianceMat;
    arma::mat Mon = Monitor;

    double Acceptance = ACC;
    int t_iter = 0;

    for (int iter = 0; iter < Iterations; iter++) {

        if ((iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                << ", LP: " << as<double>(Mo0["LP"]) << std::endl;
        }

        arma::vec prop = Bproposal(auxiliary, gr0, epsilon);

        List Mo1 = Model(wrap(prop), Data);

        double LP0 = as<double>(Mo0["LP"]);
        double LP1 = as<double>(Mo1["LP"]);

        double alpha = 1.0 / (std::exp(LP0 - LP1) + 1.0);
        double u = R::runif(0, 1);

        if (u < alpha) {
            Mo0 = Mo1;
            auxiliary = prop;
            gr0 = grad(Model, Data, prop, h);
            Acceptance += 1.0;
        }

        if ((iter + 1) % Thinning == 0) {
            t_iter = (iter + 1) / Thinning - 1;

            thinned.row(t_iter) = as<arma::vec>(Mo0["parm"]).t();
            postpred.row(t_iter) = as<arma::vec>(Mo0["yhat"]).t();
            Dev.row(t_iter) = as<arma::vec>(Mo0["Dev"]).t();
            Mon.row(t_iter) = as<arma::vec>(Mo0["Monitor"]).t();

            if (iter < Adapt) {
                arma::vec sigma = sqrt(diagvec(epsilon * epsilon.t()));
                bool accept = (u < alpha);

                epsilon = sigmaEstimate(
                    thinned.rows(0, t_iter),
                    epsilon,
                    t_iter + 1,
                    sigma,
                    accept
                );
            }
        }
    }

    return List::create(
        Named("Acceptance") = Acceptance / Iterations,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("thinned") = thinned,
        Named("postpred") = postpred
    );
}

// [[Rcpp::export]]
List ohss(Function Model, List Data, int Iterations, int Status,
          arma::vec InitialValues, int Thinning, double ACC,
          arma::mat DevianceMat, int LIV, arma::mat Monitor,
          List Mo0, arma::mat samples, arma::mat PPD, int Adapt) {

    RNGScope scope;

    arma::vec auxiliary = InitialValues;

    arma::mat thinned = samples;
    arma::mat postpred = PPD;
    arma::mat Dev = DevianceMat;
    arma::mat Mon = Monitor;

    arma::mat post(samples.n_rows, samples.n_cols, fill::zeros);

    int t_iter = 0;
    double Acceptance = ACC;

    double w = 0.05;
    double tuning = 1.0;
    double edge_scale = 5.0;

    arma::mat VarCov = eye(LIV, LIV);
    arma::mat VarCov2 = VarCov;

    arma::vec eigval;
    arma::mat eigvec;

    eig_sym(eigval, eigvec, VarCov);

    int decomp_freq = std::max((int)std::floor(Iterations / Thinning / 100.0), 10);

    for (int iter = 0; iter < Iterations; iter++) {

        if ((iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                << ", LP: " << as<double>(Mo0["LP"]) << std::endl;
        }

        // Update covariance
        if ((iter % decomp_freq == 0) && (iter > 1) && (iter < Adapt)) {
            arma::mat temp = post.rows(0, iter - 1);
            VarCov2 = (temp.t() * temp) / (iter - 1);
            eig_sym(eigval, eigvec, VarCov2);
        }

        arma::vec V_eig = eigval;
        arma::mat S_eig = eigvec;

        if (R::runif(0, 1) < w) {
            V_eig.fill(tuning);
            S_eig.eye();
        }

        double y_slice = as<double>(Mo0["LP"]) - R::rexp(1.0);

        arma::vec L = -arma::randu<arma::vec>(LIV);
        arma::vec U = L + 1.0;

        arma::vec prop = auxiliary;

        bool accepted = false;

        while (!accepted) {

            arma::vec wt = L + (U - L) % arma::randu<arma::vec>(LIV);

            arma::vec v = S_eig * (edge_scale * (wt % V_eig));
            prop = auxiliary + v;

            List Mo1 = Model(wrap(prop), Data);

            double LP1 = as<double>(Mo1["LP"]);

            if (LP1 >= y_slice) {
                auxiliary = prop;
                Mo0 = Mo1;
                accepted = true;
                break;
            }

            // shrink interval
            for (int j = 0; j < LIV; j++) {
                if (wt[j] < 0) L[j] = wt[j];
                else U[j] = wt[j];
            }

            if (all(abs(wt) < 1e-12)) {
                accepted = true;
                break;
            }
        }

        if (iter < Adapt) {
            post.row(iter) = auxiliary.t();
            VarCov = VarCov2;
        }

        Acceptance += 1.0;

        if ((iter + 1) % Thinning == 0) {
            t_iter = (iter + 1) / Thinning - 1;

            thinned.row(t_iter) = as<arma::vec>(Mo0["parm"]).t();
            postpred.row(t_iter) = as<arma::vec>(Mo0["yhat"]).t();
            Dev.row(t_iter) = as<arma::vec>(Mo0["Dev"]).t();
            Mon.row(t_iter) = as<arma::vec>(Mo0["Monitor"]).t();
        }
    }

    return List::create(
        Named("Acceptance") = Acceptance / Iterations,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("thinned") = thinned,
        Named("postpred") = postpred
    );
}

// [[Rcpp::export]]
List nuts(Function Model, List Data, int Iterations, int Status,
          arma::vec InitialValues, int Thinning, arma::mat thinned,
          arma::mat postpred, arma::mat Dev, arma::mat Mon, List Mo0,
          double h = 1e-6, int max_depth = 10, int Adapt = 1000,
          double target = 0.65) {

    RNGScope scope;

    arma::vec theta = as<arma::vec>(Mo0["parm"]);
    int dim = theta.n_elem;
    int t_iter = 0;

    double epsilon = FindReasonableEpsilon(theta, Model, Data, h);

    double mu = log(10.0 * epsilon);
    double eps_bar = epsilon;
    double Hbar = 0.0;

    double gamma = 0.05;
    double t0 = 10.0;
    double kappa = 0.75;

    double Acceptance = 0.0;

    List Mo = clone(Mo0);

    for (int iter = 0; iter < Iterations; iter++) {

        if ((iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                << ", LP: " << as<double>(Mo["LP"]) << std::endl;
        }

        arma::vec r0 = arma::randn(dim);

        double joint0 = as<double>(Mo["LP"]) - 0.5 * dot(r0, r0);
        double logu = joint0 - R::rexp(1.0);

        arma::vec theta_minus = theta;
        arma::vec theta_plus = theta;
        arma::vec r_minus = r0;
        arma::vec r_plus = r0;

        arma::vec theta_new = theta;
        arma::vec grad_current = gradC(Model, Data, theta, h);

        int n = 1;
        int s = 1;
        int j = 0;

        double alpha = 0.0;
        int n_alpha = 0;

        while (s == 1 && j < max_depth) {

            int v = (R::runif(0, 1) < 0.5) ? -1 : 1;

            List BT;

            if (v == -1) {
                BT = BuildTree(theta_minus, r_minus, grad_current,
                               logu, v, j, epsilon, theta, joint0, Model, Data, h);

                theta_minus = as<arma::vec>(BT["theta_minus"]);
                r_minus = as<arma::vec>(BT["r_minus"]);
            }
            else {
                BT = BuildTree(theta_plus, r_plus, grad_current,
                               logu, v, j, epsilon, theta, joint0, Model, Data, h);

                theta_plus = as<arma::vec>(BT["theta_plus"]);
                r_plus = as<arma::vec>(BT["r_plus"]);
            }

            arma::vec theta_prop = as<arma::vec>(BT["theta"]);
            int n_prop = as<int>(BT["n"]);
            int s_prop = as<int>(BT["s"]);

            if (s_prop == 1 && R::runif(0, 1) < (double)n_prop / n) {
                theta_new = theta_prop;
            }

            n += n_prop;

            s = s_prop &&
                (dot(theta_plus - theta_minus, r_minus) >= 0) &&
                (dot(theta_plus - theta_minus, r_plus) >= 0);

            alpha += as<double>(BT["alpha"]);
            n_alpha += as<int>(BT["n_alpha"]);

            j++;
        }

        double accept_rate = (n_alpha > 0) ? alpha / n_alpha : 0.0;
        Acceptance += accept_rate;

        if (iter < Adapt) {
            double m = iter + 1.0;

            Hbar = (1.0 - 1.0 / (m + t0)) * Hbar +
                (1.0 / (m + t0)) * (target - accept_rate);

            double log_eps = mu - (std::sqrt(m) / gamma) * Hbar;
            epsilon = std::exp(log_eps);

            double eta = std::pow(m, -kappa);
            eps_bar = std::exp(eta * log_eps + (1.0 - eta) * std::log(eps_bar));
        }
        else {
            epsilon = eps_bar;
        }

        theta = theta_new;
        Mo = Model(theta, Data);
        grad_current = gradC(Model, Data, theta, h);

        if ((iter + 1) % Thinning == 0) {
            t_iter = (iter + 1) / Thinning - 1;

            thinned.row(t_iter) = as<arma::vec>(Mo["parm"]).t();
            postpred.row(t_iter) = as<arma::vec>(Mo["yhat"]).t();
            Dev.row(t_iter) = as<arma::vec>(Mo["Dev"]).t();
            Mon.row(t_iter) = as<arma::vec>(Mo["Monitor"]).t();
        }
    }

    return List::create(
        Named("Acceptance") = Acceptance / Iterations,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("thinned") = thinned,
        Named("postpred") = postpred,
        Named("step_size") = epsilon
    );
}

// ====================================LA Algorithm====================================

arma::mat rmvnorm_arma(int n, const arma::vec& mu, const arma::mat& Sigma) {
  int d = mu.n_elem;
  arma::mat Z = arma::randn(n, d);
  arma::mat C = arma::chol(Sigma);
  return Z * C + arma::repmat(mu.t(), n, 1);
}

arma::vec dmvnorm(const arma::mat& X, const arma::vec& mean, const arma::mat& Sigma) {

    int n = X.n_rows;
    int d = X.n_cols;

    arma::mat Xc = X.each_row() - mean.t();

    // FIX: use Cholesky for both log-det and quadratic form -- numerically stable
    arma::mat C = arma::chol(Sigma);
    double logdetSigma = 2.0 * arma::sum(arma::log(arma::diagvec(C)));

    // Solve C^T C x = Xc^T via forward/back substitution
    arma::mat solved = arma::solve(arma::trimatl(C.t()), Xc.t());

    arma::vec out(n);
    for (int i = 0; i < n; i++) {
        double quad = arma::dot(solved.col(i), solved.col(i));
        out[i] = -0.5 * (d * std::log(2.0 * M_PI) + logdetSigma + quad);
    }

    return out;
}

List psis_smooth_weights_cpp(arma::vec log_w, double tail_frac = 0.2) {

    int S = log_w.n_elem;

    // Stabilize
    double maxlw = log_w.max();
    arma::vec lw = log_w - maxlw;

    // Convert to weights
    arma::vec w = exp(lw);

    // Sort descending
    arma::uvec ord = sort_index(w, "descend");
    arma::vec w_sorted = w.elem(ord);

    // Tail size
    int M = std::max(5, (int)std::floor(tail_frac * S));
    arma::vec w_tail = w_sorted.subvec(0, M - 1);

    double threshold = w_tail[M - 1];
    arma::vec excess = w_tail - threshold;

    // ===== Hill estimator =====
    double k_hat = 0.0;

    if (all(excess <= 0)) {
        k_hat = 0.0;
    }
    else {
        arma::vec ratio = w_tail / threshold;
        k_hat = mean(log(ratio + 1e-12));
    }

    // ===== Smooth tail =====
    arma::vec ranks = linspace<vec>(1, M, M);
    arma::vec p = (ranks - 0.5) / M;

    arma::vec smoothed(M);

    if (std::abs(k_hat) < 1e-8) {
        smoothed = threshold - log(1 - p);
    }
    else {
        smoothed = threshold * pow(1 - p, -k_hat);
    }

    // Replace tail
    w_sorted.subvec(0, M - 1) = smoothed;

    // Restore original order
    arma::vec w_new(S);
    w_new.elem(ord) = w_sorted;

    // Normalize
    w_new = w_new / sum(w_new);

    arma::vec log_w_new = log(w_new + 1e-300);

    return List::create(
        Named("w") = w_new,
        Named("log_w") = log_w_new,
        Named("pareto_k") = k_hat
    );
}

// [[Rcpp::export]]
List samplingImportanceResampling(arma::vec MAP, arma::mat VarCov, Function Model,
                                  List Data, int iterations) {

    RNGScope scope;

    int d = MAP.n_elem;

    // ===== Sample from proposal =====
    arma::mat theta = rmvnorm_arma(iterations, MAP, VarCov);
    arma::mat posterior(iterations, d);

    // ===== Initial evaluation =====
    List tmp0 = Model(wrap(MAP), Data);

    int ydim = as<arma::vec>(tmp0["yhat"]).n_elem;
    int mdim = as<arma::vec>(tmp0["Monitor"]).n_elem;

    arma::mat yhat(iterations, ydim);
    arma::mat Monitor(iterations, mdim);
    arma::vec LP(iterations);
    arma::vec Dev(iterations);

    // ===== Evaluate posterior =====
    for (int i = 0; i < iterations; i++) {
        List tmp = Model(wrap(theta.row(i).t()), Data);

        posterior.row(i) = as<arma::vec>(tmp["parm"]).t();
        LP[i] = as<double>(tmp["LP"]);
        Dev[i] = as<double>(tmp["Dev"]);
        yhat.row(i) = as<arma::vec>(tmp["yhat"]).t();
        Monitor.row(i) = as<arma::vec>(tmp["Monitor"]).t();
    }

    // ===== Proposal density =====
    arma::vec LMVN = dmvnorm(theta, MAP, VarCov);

    // ===== Raw log-weights =====
    arma::vec log_w = LP - LMVN;
    log_w = log_w - log_w.max();

    // Handle bad values
    arma::uvec bad = find_nonfinite(log_w);
    if (bad.n_elem > 0) {
        double min_good = log_w.elem(find_finite(log_w)).min();
        log_w.elem(bad).fill(min_good);
    }

    arma::vec w_raw = exp(log_w);
    w_raw /= sum(w_raw);

    // ===== PSIS =====
    List psis = psis_smooth_weights_cpp(log_w);

    arma::vec w_psis = as<arma::vec>(psis["w"]);
    arma::vec log_w_psis = as<arma::vec>(psis["log_w"]);
    double pareto_k = as<double>(psis["pareto_k"]);

    // ===== Resampling =====
    arma::vec u = arma::randu(iterations);
    arma::vec cdf = cumsum(w_psis);

    arma::uvec indices(iterations);

    for (int i = 0; i < iterations; i++) {
        indices[i] = arma::as_scalar(find(cdf >= u[i], 1, "first"));
    }

    arma::mat posterior_rs = posterior.rows(indices);
    arma::mat yhat_rs = yhat.rows(indices);
    arma::mat Mon_rs = Monitor.rows(indices);
    arma::vec LP_rs = LP.elem(indices);
    arma::vec Dev_rs = Dev.elem(indices);

    // ===== Output =====
    return List::create(Named("indices") = indices + 1,
                        Named("posterior") = posterior_rs,
                        Named("yhat") = yhat_rs,
                        Named("Monitor") = Mon_rs,
                        Named("LP") = LP_rs,
                        Named("Dev") = Dev_rs,
                        Named("weights_raw") = w_raw,
                        Named("log_w_raw") = log_w,
                        Named("weights_psis") = w_psis,
                        Named("log_w_psis") = log_w_psis,
                        Named("pareto_k") = pareto_k
    );
}

// ====================================VB Algorithms====================================

// The variational algorithms below work with the unnormalised log posterior returned
// in Model(...)["LP"].  Numerical derivatives are scale-aware and all covariance /
// precision matrices are projected back to the positive-definite cone when needed.

inline double vb_lp(Function Model, List Data, const arma::vec& x) {
    List out = Model(wrap(x), Data);
    double value = as<double>(out["LP"]);
    return std::isfinite(value) ? value : -std::numeric_limits<double>::infinity();
}

arma::vec vb_steps(const arma::vec& x, double h) {
    arma::vec out = h * arma::max(arma::abs(x), arma::ones<arma::vec>(x.n_elem));
    out.transform([](double v) { return std::max(v, 1e-8); });
    return out;
}

arma::vec vb_grad(Function Model, List Data, const arma::vec& x, double h = 1e-4) {
    int d = x.n_elem;
    arma::vec out(d, fill::zeros);
    arma::vec hs = vb_steps(x, h);

    for (int j = 0; j < d; ++j) {
        arma::vec xp = x;
        arma::vec xm = x;
        xp[j] += hs[j];
        xm[j] -= hs[j];
        double fp = vb_lp(Model, Data, xp);
        double fm = vb_lp(Model, Data, xm);

        if (std::isfinite(fp) && std::isfinite(fm)) {
            out[j] = (fp - fm) / (2.0 * hs[j]);
        } else {
            // One-sided fallback around the current point.
            double f0 = vb_lp(Model, Data, x);
            if (std::isfinite(fp) && std::isfinite(f0)) {
                out[j] = (fp - f0) / hs[j];
            } else if (std::isfinite(fm) && std::isfinite(f0)) {
                out[j] = (f0 - fm) / hs[j];
            }
        }
        if (!std::isfinite(out[j])) out[j] = 0.0;
    }
    return out;
}

arma::mat vb_hessian(Function Model, List Data, const arma::vec& x,
                     double h = 1e-4) {
    int d = x.n_elem;
    arma::mat H(d, d, fill::zeros);
    arma::vec hs = vb_steps(x, h);
    double f0 = vb_lp(Model, Data, x);

    for (int i = 0; i < d; ++i) {
        arma::vec xp = x;
        arma::vec xm = x;
        xp[i] += hs[i];
        xm[i] -= hs[i];
        double fp = vb_lp(Model, Data, xp);
        double fm = vb_lp(Model, Data, xm);
        if (std::isfinite(fp) && std::isfinite(fm) && std::isfinite(f0)) {
            H(i, i) = (fp - 2.0 * f0 + fm) / (hs[i] * hs[i]);
        }

        for (int j = i + 1; j < d; ++j) {
            arma::vec xpp = x;
            arma::vec xpm = x;
            arma::vec xmp = x;
            arma::vec xmm = x;
            xpp[i] += hs[i]; xpp[j] += hs[j];
            xpm[i] += hs[i]; xpm[j] -= hs[j];
            xmp[i] -= hs[i]; xmp[j] += hs[j];
            xmm[i] -= hs[i]; xmm[j] -= hs[j];
            double fpp = vb_lp(Model, Data, xpp);
            double fpm = vb_lp(Model, Data, xpm);
            double fmp = vb_lp(Model, Data, xmp);
            double fmm = vb_lp(Model, Data, xmm);
            if (std::isfinite(fpp) && std::isfinite(fpm) &&
                std::isfinite(fmp) && std::isfinite(fmm)) {
                double hij = (fpp - fpm - fmp + fmm) /
                    (4.0 * hs[i] * hs[j]);
                H(i, j) = hij;
                H(j, i) = hij;
            }
        }
    }

    H = 0.5 * (H + H.t());
    return H;
}

arma::mat vb_make_pd(const arma::mat& A, double eig_floor = 1e-8,
                     double eig_ceiling = 1e8) {
    arma::mat S = 0.5 * (A + A.t());
    arma::vec eigval;
    arma::mat eigvec;
    bool ok = arma::eig_sym(eigval, eigvec, S);
    if (!ok || eigval.n_elem == 0 || !eigval.is_finite()) {
        return arma::eye(A.n_rows, A.n_cols);
    }
    eigval.transform([&](double v) {
        if (!std::isfinite(v)) return eig_floor;
        return std::min(std::max(v, eig_floor), eig_ceiling);
    });
    return eigvec * arma::diagmat(eigval) * eigvec.t();
}

arma::mat vb_safe_inverse(const arma::mat& A, double eig_floor = 1e-8) {
    arma::mat P = vb_make_pd(A, eig_floor);
    arma::mat out;
    bool ok = arma::inv_sympd(out, P);
    if (!ok || !out.is_finite()) {
        out = arma::pinv(P);
    }
    return vb_make_pd(out, eig_floor);
}

arma::mat vb_rmvnorm(int n, const arma::vec& mu, const arma::mat& Sigma,
                     double eig_floor = 1e-8) {
    arma::mat S = vb_make_pd(Sigma, eig_floor);
    arma::mat C;
    bool ok = arma::chol(C, S, "lower");
    if (!ok) {
        S += 1e-6 * arma::eye(S.n_rows, S.n_cols);
        arma::chol(C, S, "lower");
    }
    arma::mat Z = arma::randn(n, mu.n_elem);
    return Z * C.t() + arma::repmat(mu.t(), n, 1);
}

inline double vb_log_q_gaussian(const arma::vec& x, const arma::vec& mu,
                                const arma::mat& Sigma) {
    arma::mat S = vb_make_pd(Sigma, 1e-10);
    arma::mat C;
    bool ok = arma::chol(C, S, "lower");
    if (!ok) return -std::numeric_limits<double>::infinity();
    arma::vec z = arma::solve(arma::trimatl(C), x - mu);
    double logdet = 2.0 * arma::sum(arma::log(C.diag()));
    return -0.5 * (x.n_elem * std::log(2.0 * M_PI) + logdet + arma::dot(z, z));
}

inline double vb_log_diag_normal(const arma::vec& x, const arma::vec& mu,
                                 const arma::vec& log_sd) {
    arma::vec z = (x - mu) % arma::exp(-log_sd);
    return -0.5 * x.n_elem * std::log(2.0 * M_PI) -
        arma::sum(log_sd) - 0.5 * arma::dot(z, z);
}

void vb_store_draw(Function Model, List Data, const arma::vec& x, int row,
                   arma::mat& thinned, arma::mat& postpred,
                   arma::mat& Dev, arma::mat& Mon,
                   arma::mat& latent, arma::vec& LP) {
    List out = Model(wrap(x), Data);
    thinned.row(row) = as<arma::vec>(out["parm"]).t();
    postpred.row(row) = as<arma::vec>(out["yhat"]).t();
    Dev.row(row) = as<arma::vec>(out["Dev"]).t();
    Mon.row(row) = as<arma::vec>(out["Monitor"]).t();
    latent.row(row) = x.t();
    LP[row] = as<double>(out["LP"]);
}

bool vb_convergence_update(double criterion, double tolerance, int iter,
                           int min_iter, int patience, int& stable_count) {
    if (tolerance <= 0.0 || iter + 1 < min_iter || !std::isfinite(criterion)) {
        stable_count = 0;
        return false;
    }
    if (criterion < tolerance) ++stable_count;
    else stable_count = 0;
    return stable_count >= patience;
}

arma::mat vb_reconstruct_lbfgs_hessian(
    const std::vector<arma::vec>& s_hist,
    const std::vector<arma::vec>& y_hist,
    int d, double damping, double eig_floor) {

    if (s_hist.empty()) return arma::eye(d, d);

    const arma::vec& sl = s_hist.back();
    const arma::vec& yl = y_hist.back();
    double ys_last = arma::dot(yl, sl);
    double gamma = 1.0;
    if (ys_last > 1e-12) {
        gamma = arma::dot(yl, yl) / ys_last;
    }
    gamma = std::min(std::max(gamma, 1e-6), 1e6);
    arma::mat B = gamma * arma::eye(d, d);

    for (std::size_t k = 0; k < s_hist.size(); ++k) {
        arma::vec s = s_hist[k];
        arma::vec y = y_hist[k];
        arma::vec Bs = B * s;
        double sBs = arma::dot(s, Bs);
        double ys = arma::dot(y, s);
        if (!(sBs > 1e-12) || !std::isfinite(sBs) || !std::isfinite(ys)) continue;

        if (ys < damping * sBs) {
            double denom = sBs - ys;
            if (std::abs(denom) > 1e-12) {
                double theta = (1.0 - damping) * sBs / denom;
                theta = std::min(std::max(theta, 0.0), 1.0);
                y = theta * y + (1.0 - theta) * Bs;
                ys = arma::dot(y, s);
            }
        }
        if (!(ys > 1e-12)) continue;

        B -= (Bs * Bs.t()) / sBs;
        B += (y * y.t()) / ys;
        B = vb_make_pd(B, eig_floor);
    }
    return B;
}

// [[Rcpp::export]]
List sagva(Function Model, List Data, int Iterations, int Status,
           arma::vec InitialValues, arma::mat InitialCov, int Thinning,
           arma::mat thinned, arma::mat postpred,
           arma::mat Dev, arma::mat Mon,
           double h = 1e-4,
           double learning_rate = -1.0,
           double Stop_Tolerance = 1e-4,
           int Min_Iterations = 200,
           int Patience = 50,
           double eig_floor = 1e-8) {

    RNGScope scope;
    int d = InitialValues.n_elem;
    int n_draws = thinned.n_rows;

    arma::vec m = InitialValues;
    arma::mat V = vb_make_pd(InitialCov, eig_floor);
    arma::mat P = vb_safe_inverse(V, eig_floor);
    arma::vec a(d, fill::zeros);
    arma::vec z = m;

    arma::vec avg_g(d, fill::zeros);
    arma::vec avg_x(d, fill::zeros);
    arma::mat avg_precision(d, d, fill::zeros);
    int avg_count = 0;

    //double w = learning_rate > 0.0 ? learning_rate : 1.0 / std::sqrt((double)Iterations);
    //w = std::min(std::max(w, 1e-5), 1.0);
    double w = 0.0;
    double base_w = learning_rate > 0.0 ? learning_rate : 0.10;

    double criterion = std::numeric_limits<double>::infinity();
    double criterion_ema = criterion;
    int stable_count = 0;
    int completed = 0;
    int invalid_updates = 0;
    bool converged = false;

    for (int iter = 0; iter < Iterations; ++iter) {
        w = base_w * std::pow(10.0 / (10.0 + iter), 0.60);
        w = std::min(std::max(w, 1e-5), 1.0);

        arma::vec old_m = m;
        arma::mat old_V = V;

        arma::vec x = vb_rmvnorm(1, m, V, eig_floor).row(0).t();
        arma::vec g = vb_grad(Model, Data, x, h);
        arma::mat H = vb_hessian(Model, Data, x, h);

        if (!g.is_finite() || !H.is_finite()) {
            ++invalid_updates;
            V *= 0.5;
            V = vb_make_pd(V, eig_floor);
            P = vb_safe_inverse(V, eig_floor);
            continue;
        }

        a = (1.0 - w) * a + w * g;
        z = (1.0 - w) * z + w * x;
        P = (1.0 - w) * P - w * H;
        P = vb_make_pd(P, eig_floor);
        V = vb_safe_inverse(P, eig_floor);
        m = V * a + z;

        if (!m.is_finite() || !V.is_finite()) {
            ++invalid_updates;
            m = old_m;
            V = 0.5 * old_V;
            V = vb_make_pd(V, eig_floor);
            P = vb_safe_inverse(V, eig_floor);
            continue;
        }

        if (iter >= Iterations / 2) {
            ++avg_count;
            avg_g += (g - avg_g) / (double)avg_count;
            avg_x += (x - avg_x) / (double)avg_count;
            arma::mat negH = -H;
            avg_precision += (negH - avg_precision) / (double)avg_count;
        }

        double dm = arma::norm(m - old_m, 2) / (1.0 + arma::norm(old_m, 2));
        double dV = arma::norm(V - old_V, "fro") / (1.0 + arma::norm(old_V, "fro"));
        criterion = std::max(dm, dV);
        if (!std::isfinite(criterion_ema)) criterion_ema = criterion;
        else criterion_ema = 0.95 * criterion_ema + 0.05 * criterion;

        completed = iter + 1;
        if (Status > 0 && (iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                  << ", LP(mean): " << vb_lp(Model, Data, m)
                  << ", convergence criterion: " << criterion_ema << std::endl;
        }

        if (vb_convergence_update(criterion_ema, Stop_Tolerance, iter,
                                  Min_Iterations, Patience, stable_count)) {
            converged = true;
            break;
        }
    }

    arma::vec final_m = m;
    arma::mat final_V = V;
    if (avg_count > 1) {
        arma::mat final_P = vb_make_pd(avg_precision, eig_floor);
        final_V = vb_safe_inverse(final_P, eig_floor);
        final_m = final_V * avg_g + avg_x;
        if (!final_m.is_finite()) final_m = m;
    }

    arma::mat latent(n_draws, d, fill::zeros);
    arma::vec LP(n_draws, fill::zeros);
    arma::vec log_q(n_draws, fill::zeros);
    arma::mat draws = vb_rmvnorm(n_draws, final_m, final_V, eig_floor);
    for (int i = 0; i < n_draws; ++i) {
        arma::vec x = draws.row(i).t();
        vb_store_draw(Model, Data, x, i, thinned, postpred, Dev, Mon, latent, LP);
        log_q[i] = vb_log_q_gaussian(x, final_m, final_V);
    }
    arma::vec elbo_samples = LP - log_q;

    return List::create(
        Named("thinned") = thinned,
        Named("postpred") = postpred,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("latent") = latent,
        Named("LP") = LP,
        Named("log_q") = log_q,
        Named("ELBO_samples") = elbo_samples,
        Named("LowerBound") = arma::mean(elbo_samples),
        Named("mean") = final_m,
        Named("VarCov") = final_V,
        Named("criterion") = criterion_ema,
        Named("n_iter") = completed,
        Named("converged") = converged,
        Named("invalid_updates") = invalid_updates,
        Named("learning_rate") = w
    );
}

// [[Rcpp::export]]
List qnsagva(Function Model, List Data, int Iterations, int Status,
             arma::vec InitialValues, arma::mat InitialCov, int Thinning,
             arma::mat thinned, arma::mat postpred,
             arma::mat Dev, arma::mat Mon,
             int memory = 10,
             double h = 1e-4,
             double damping = 0.2,
             double learning_rate = -1.0,
             double Stop_Tolerance = 1e-4,
             int Min_Iterations = 200,
             int Patience = 50,
             double eig_floor = 1e-8) {

    RNGScope scope;
    int d = InitialValues.n_elem;
    int n_draws = thinned.n_rows;
    memory = std::max(memory, 1);

    arma::vec m = InitialValues;
    arma::mat V = vb_make_pd(InitialCov, eig_floor);
    arma::mat P = vb_safe_inverse(V, eig_floor);
    arma::vec a(d, fill::zeros);
    arma::vec z = m;

    std::vector<arma::vec> s_hist;
    std::vector<arma::vec> y_hist;
    arma::vec previous_x;
    arma::vec previous_g;
    bool have_previous = false;

    arma::vec avg_g(d, fill::zeros);
    arma::vec avg_x(d, fill::zeros);
    arma::mat avg_precision(d, d, fill::zeros);
    int avg_count = 0;

    //double w = learning_rate > 0.0 ? learning_rate : 1.0 / std::sqrt((double)Iterations);
    //w = std::min(std::max(w, 1e-5), 1.0);
    double base_w = learning_rate > 0.0 ? learning_rate : 0.10;
    double w = 0.0;

    double criterion = std::numeric_limits<double>::infinity();
    double criterion_ema = criterion;
    int stable_count = 0;
    int completed = 0;
    int invalid_updates = 0;
    int accepted_pairs = 0;
    bool converged = false;

    for (int iter = 0; iter < Iterations; ++iter) {
        double w = base_w * std::pow(10.0 / (10.0 + iter), 0.60);
        w = std::min(std::max(w, 1e-5), 1.0);

        arma::vec old_m = m;
        arma::mat old_V = V;

        arma::vec x = vb_rmvnorm(1, m, V, eig_floor).row(0).t();
        arma::vec g = vb_grad(Model, Data, x, h);
        if (!g.is_finite()) {
            ++invalid_updates;
            V *= 0.5;
            V = vb_make_pd(V, eig_floor);
            P = vb_safe_inverse(V, eig_floor);
            continue;
        }

        if (have_previous) {
            arma::vec s = x - previous_x;
            arma::vec y = previous_g - g; // Hessian of -log p times s near a mode.
            double ys = arma::dot(y, s);
            double scale = arma::norm(s, 2) * arma::norm(y, 2);
            if (ys > 1e-10 * std::max(1.0, scale)) {
                s_hist.push_back(s);
                y_hist.push_back(y);
                if ((int)s_hist.size() > memory) {
                    s_hist.erase(s_hist.begin());
                    y_hist.erase(y_hist.begin());
                }
                ++accepted_pairs;
            }
        }
        previous_x = x;
        previous_g = g;
        have_previous = true;

        arma::mat B = vb_reconstruct_lbfgs_hessian(
            s_hist, y_hist, d, damping, eig_floor);
        if (s_hist.empty()) B = P;

        a = (1.0 - w) * a + w * g;
        z = (1.0 - w) * z + w * x;
        P = (1.0 - w) * P + w * B;
        P = vb_make_pd(P, eig_floor);
        V = vb_safe_inverse(P, eig_floor);
        m = V * a + z;

        if (!m.is_finite() || !V.is_finite()) {
            ++invalid_updates;
            m = old_m;
            V = 0.5 * old_V;
            V = vb_make_pd(V, eig_floor);
            P = vb_safe_inverse(V, eig_floor);
            continue;
        }

        if (iter >= Iterations / 2) {
            ++avg_count;
            avg_g += (g - avg_g) / (double)avg_count;
            avg_x += (x - avg_x) / (double)avg_count;
            avg_precision += (B - avg_precision) / (double)avg_count;
        }

        double dm = arma::norm(m - old_m, 2) / (1.0 + arma::norm(old_m, 2));
        double dV = arma::norm(V - old_V, "fro") / (1.0 + arma::norm(old_V, "fro"));
        criterion = std::max(dm, dV);
        if (!std::isfinite(criterion_ema)) criterion_ema = criterion;
        else criterion_ema = 0.95 * criterion_ema + 0.05 * criterion;

        completed = iter + 1;
        if (Status > 0 && (iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                  << ", LP(mean): " << vb_lp(Model, Data, m)
                  << ", curvature pairs: " << s_hist.size()
                  << ", convergence criterion: " << criterion_ema << std::endl;
        }

        if (vb_convergence_update(criterion_ema, Stop_Tolerance, iter,
                                  Min_Iterations, Patience, stable_count)) {
            converged = true;
            break;
        }
    }

    arma::vec final_m = m;
    arma::mat final_V = V;
    if (avg_count > 1) {
        arma::mat final_P = vb_make_pd(avg_precision, eig_floor);
        final_V = vb_safe_inverse(final_P, eig_floor);
        final_m = final_V * avg_g + avg_x;
        if (!final_m.is_finite()) final_m = m;
    }

    arma::mat latent(n_draws, d, fill::zeros);
    arma::vec LP(n_draws, fill::zeros);
    arma::vec log_q(n_draws, fill::zeros);
    arma::mat draws = vb_rmvnorm(n_draws, final_m, final_V, eig_floor);
    for (int i = 0; i < n_draws; ++i) {
        arma::vec x = draws.row(i).t();
        vb_store_draw(Model, Data, x, i, thinned, postpred, Dev, Mon, latent, LP);
        log_q[i] = vb_log_q_gaussian(x, final_m, final_V);
    }
    arma::vec elbo_samples = LP - log_q;

    return List::create(
        Named("thinned") = thinned,
        Named("postpred") = postpred,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("latent") = latent,
        Named("LP") = LP,
        Named("log_q") = log_q,
        Named("ELBO_samples") = elbo_samples,
        Named("LowerBound") = arma::mean(elbo_samples),
        Named("mean") = final_m,
        Named("VarCov") = final_V,
        Named("criterion") = criterion_ema,
        Named("n_iter") = completed,
        Named("converged") = converged,
        Named("invalid_updates") = invalid_updates,
        Named("accepted_curvature_pairs") = accepted_pairs,
        Named("learning_rate") = w
    );
}

// Run one Hamiltonian variational trajectory.  The forward momentum is standard
// normal; the learned reverse model is diagonal Gaussian with mean
// a_r * grad log p(z_t) + c_r * z_t + b_r.
double vb_hvi_bound(Function Model, List Data,
                    const arma::vec& theta,
                    const std::vector<arma::vec>& eps0,
                    const std::vector<std::vector<arma::vec> >& momentum,
                    const std::vector<arma::vec>& a_r,
                    const std::vector<arma::vec>& c_r,
                    const std::vector<arma::vec>& b_r,
                    const std::vector<arma::vec>& log_sd_r,
                    int T_mcmc, int L_leapfrog, double h,
                    int& invalid_count) {

    int d = theta.n_elem / 2;
    arma::vec mu = theta.subvec(0, d - 1);
    arma::vec log_sd = theta.subvec(d, 2 * d - 1);
    double epsilon = std::exp(theta[2 * d]);
    double total = 0.0;

    for (std::size_t k = 0; k < eps0.size(); ++k) {
        arma::vec z = mu + arma::exp(log_sd) % eps0[k];
        double lp_prev = vb_lp(Model, Data, z);
        if (!std::isfinite(lp_prev)) {
            ++invalid_count;
            total += -1e12;
            continue;
        }
        double Lk = lp_prev - vb_log_diag_normal(z, mu, log_sd);

        for (int t = 0; t < T_mcmc; ++t) {
            arma::vec v0 = momentum[k][t];
            double log_q_v0 = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * arma::dot(v0, v0);
            arma::vec z_new = z;
            arma::vec v_new = v0;
            arma::vec g = vb_grad(Model, Data, z_new, h);

            bool valid = g.is_finite();
            for (int l = 0; l < L_leapfrog && valid; ++l) {
                v_new += 0.5 * epsilon * g;
                z_new += epsilon * v_new;
                g = vb_grad(Model, Data, z_new, h);
                if (!g.is_finite() || !z_new.is_finite() || !v_new.is_finite()) valid = false;
                if (valid) v_new += 0.5 * epsilon * g;
            }

            double lp_new = valid ? vb_lp(Model, Data, z_new) :
                -std::numeric_limits<double>::infinity();
            if (!std::isfinite(lp_new)) {
                ++invalid_count;
                Lk += -1e12;
                break;
            }

            arma::vec mu_r = a_r[t] % g + c_r[t] % z_new + b_r[t];
            double log_r = vb_log_diag_normal(v_new, mu_r, log_sd_r[t]);
            Lk += lp_new + log_r - lp_prev - log_q_v0;
            z = z_new;
            lp_prev = lp_new;
        }
        total += Lk;
    }
    return total / (double)eps0.size();
}

// [[Rcpp::export]]
List mcvi(Function Model, List Data, int Iterations, int Status,
          arma::mat thinned, arma::mat postpred,
          arma::mat Dev, arma::mat Mon,
          arma::vec InitialValues, arma::mat InitialCov,
          int Thinning,
          int T_mcmc = 3,
          int L_leapfrog = 3,
          int K_rb = 5,
          double h = 1e-4,
          double learning_rate = 5e-3,
          double spsa_scale = 1e-3,
          double initial_epsilon = -1.0,
          double Stop_Tolerance = 1e-4,
          int Min_Iterations = 200,
          int Patience = 50,
          double grad_clip = 100.0) {

    RNGScope scope;
    int d = InitialValues.n_elem;
    int n_draws = thinned.n_rows;
    T_mcmc = std::max(T_mcmc, 1);
    L_leapfrog = std::max(L_leapfrog, 1);
    K_rb = std::max(K_rb, 1);

    arma::vec mu_q = InitialValues;
    arma::vec init_var = vb_make_pd(InitialCov, 1e-8).diag();
    init_var.transform([](double v) { return std::max(v, 1e-10); });
    arma::vec log_sd_q = 0.5 * arma::log(init_var);

    double epsilon = initial_epsilon;
    if (!(epsilon > 0.0) || !std::isfinite(epsilon)) {
        epsilon = FindReasonableEpsilon(InitialValues, Model, Data, h);
    }
    epsilon = std::min(std::max(epsilon, 1e-5), 1.0);

    arma::vec theta(2 * d + 1, fill::zeros);
    theta.subvec(0, d - 1) = mu_q;
    theta.subvec(d, 2 * d - 1) = log_sd_q;
    theta[2 * d] = std::log(epsilon);

    std::vector<arma::vec> a_r(T_mcmc, arma::zeros(d));
    std::vector<arma::vec> c_r(T_mcmc, arma::zeros(d));
    std::vector<arma::vec> b_r(T_mcmc, arma::zeros(d));
    std::vector<arma::vec> log_sd_r(T_mcmc, arma::zeros(d));

    arma::vec adam_m(theta.n_elem, fill::zeros);
    arma::vec adam_v(theta.n_elem, fill::zeros);
    std::vector<arma::vec> m_ar(T_mcmc, arma::zeros(d)), v_ar(T_mcmc, arma::zeros(d));
    std::vector<arma::vec> m_cr(T_mcmc, arma::zeros(d)), v_cr(T_mcmc, arma::zeros(d));
    std::vector<arma::vec> m_br(T_mcmc, arma::zeros(d)), v_br(T_mcmc, arma::zeros(d));
    std::vector<arma::vec> m_lr(T_mcmc, arma::zeros(d)), v_lr(T_mcmc, arma::zeros(d));

    double beta1 = 0.9, beta2 = 0.999, adam_eps = 1e-8;
    auto adam_update = [&](arma::vec& par, const arma::vec& gradv,
                           arma::vec& m1, arma::vec& m2, int t, double lr) {
        arma::vec g = gradv;
        double ng = arma::norm(g, 2);
        if (ng > grad_clip && ng > 0.0) g *= grad_clip / ng;
        m1 = beta1 * m1 + (1.0 - beta1) * g;
        m2 = beta2 * m2 + (1.0 - beta2) * (g % g);
        arma::vec mhat = m1 / (1.0 - std::pow(beta1, t));
        arma::vec vhat = m2 / (1.0 - std::pow(beta2, t));
        par += lr * mhat / (arma::sqrt(vhat) + adam_eps);
    };

    arma::vec elbo_history(Iterations);
    elbo_history.fill(arma::datum::nan);
    double elbo_ema = -std::numeric_limits<double>::infinity();
    double criterion = std::numeric_limits<double>::infinity();
    int stable_count = 0;
    int invalid_total = 0;
    int completed = 0;
    bool converged = false;

    for (int iter = 0; iter < Iterations; ++iter) {
        double lr_now = learning_rate / std::sqrt(1.0 + iter / 100.0);

        std::vector<arma::vec> eps0(K_rb);
        std::vector<std::vector<arma::vec> > momentum(
            K_rb, std::vector<arma::vec>(T_mcmc));
        for (int k = 0; k < K_rb; ++k) {
            eps0[k] = arma::randn(d);
            for (int t = 0; t < T_mcmc; ++t) momentum[k][t] = arma::randn(d);
        }

        arma::vec delta(theta.n_elem);
        for (arma::uword j = 0; j < delta.n_elem; ++j) {
            delta[j] = R::runif(0.0, 1.0) < 0.5 ? -1.0 : 1.0;
        }
        double c_t = spsa_scale / std::pow(iter + 1.0, 0.101);
        c_t = std::max(c_t, 1e-6);
        arma::vec theta_plus = theta + c_t * delta;
        arma::vec theta_minus = theta - c_t * delta;

        // Keep scales and HMC step size in a safe numerical range.
        theta_plus.subvec(d, 2 * d - 1) =
            arma::clamp(theta_plus.subvec(d, 2 * d - 1), -10.0, 5.0);
        theta_minus.subvec(d, 2 * d - 1) =
            arma::clamp(theta_minus.subvec(d, 2 * d - 1), -10.0, 5.0);
        theta_plus[2 * d] = std::min(std::max(theta_plus[2 * d], std::log(1e-5)), std::log(1.0));
        theta_minus[2 * d] = std::min(std::max(theta_minus[2 * d], std::log(1e-5)), std::log(1.0));

        int invalid_plus = 0, invalid_minus = 0;
        double Lplus = vb_hvi_bound(Model, Data, theta_plus, eps0, momentum,
                                    a_r, c_r, b_r, log_sd_r,
                                    T_mcmc, L_leapfrog, h, invalid_plus);
        double Lminus = vb_hvi_bound(Model, Data, theta_minus, eps0, momentum,
                                     a_r, c_r, b_r, log_sd_r,
                                     T_mcmc, L_leapfrog, h, invalid_minus);
        arma::vec g_theta = ((Lplus - Lminus) / (2.0 * c_t)) * delta;

        arma::vec old_theta = theta;
        adam_update(theta, g_theta, adam_m, adam_v, iter + 1, lr_now);
        theta.subvec(d, 2 * d - 1) = arma::clamp(theta.subvec(d, 2 * d - 1), -10.0, 5.0);
        theta[2 * d] = std::min(std::max(theta[2 * d], std::log(1e-5)), std::log(1.0));

        // Current trajectories provide unbiased gradients for the reverse-model
        // parameters because those parameters do not affect the forward chain.
        arma::vec mu = theta.subvec(0, d - 1);
        arma::vec log_sd = theta.subvec(d, 2 * d - 1);
        epsilon = std::exp(theta[2 * d]);
        std::vector<arma::vec> g_ar(T_mcmc, arma::zeros(d));
        std::vector<arma::vec> g_cr(T_mcmc, arma::zeros(d));
        std::vector<arma::vec> g_br(T_mcmc, arma::zeros(d));
        std::vector<arma::vec> g_lr(T_mcmc, arma::zeros(d));
        double Lcur = 0.0;
        int invalid_cur = 0;

        for (int k = 0; k < K_rb; ++k) {
            arma::vec z = mu + arma::exp(log_sd) % eps0[k];
            double lp_prev = vb_lp(Model, Data, z);
            if (!std::isfinite(lp_prev)) {
                ++invalid_cur;
                Lcur += -1e12;
                continue;
            }
            double Lk = lp_prev - vb_log_diag_normal(z, mu, log_sd);

            for (int t = 0; t < T_mcmc; ++t) {
                arma::vec v0 = momentum[k][t];
                double log_q_v0 = -0.5 * d * std::log(2.0 * M_PI) - 0.5 * arma::dot(v0, v0);
                arma::vec z_new = z;
                arma::vec v_new = v0;
                arma::vec g = vb_grad(Model, Data, z_new, h);
                bool valid = g.is_finite();
                for (int l = 0; l < L_leapfrog && valid; ++l) {
                    v_new += 0.5 * epsilon * g;
                    z_new += epsilon * v_new;
                    g = vb_grad(Model, Data, z_new, h);
                    if (!g.is_finite() || !z_new.is_finite() || !v_new.is_finite()) valid = false;
                    if (valid) v_new += 0.5 * epsilon * g;
                }
                double lp_new = valid ? vb_lp(Model, Data, z_new) :
                    -std::numeric_limits<double>::infinity();
                if (!std::isfinite(lp_new)) {
                    ++invalid_cur;
                    Lk += -1e12;
                    break;
                }

                arma::vec mu_r = a_r[t] % g + c_r[t] % z_new + b_r[t];
                arma::vec inv_var = arma::exp(-2.0 * log_sd_r[t]);
                arma::vec resid = v_new - mu_r;
                double log_r = vb_log_diag_normal(v_new, mu_r, log_sd_r[t]);
                Lk += lp_new + log_r - lp_prev - log_q_v0;

                arma::vec score_mu = resid % inv_var;
                g_ar[t] += score_mu % g;
                g_cr[t] += score_mu % z_new;
                g_br[t] += score_mu;
                g_lr[t] += -arma::ones<arma::vec>(d) + (resid % resid) % inv_var;

                z = z_new;
                lp_prev = lp_new;
            }
            Lcur += Lk;
        }

        Lcur /= (double)K_rb;
        for (int t = 0; t < T_mcmc; ++t) {
            g_ar[t] /= (double)K_rb;
            g_cr[t] /= (double)K_rb;
            g_br[t] /= (double)K_rb;
            g_lr[t] /= (double)K_rb;
            adam_update(a_r[t], g_ar[t], m_ar[t], v_ar[t], iter + 1, lr_now);
            adam_update(c_r[t], g_cr[t], m_cr[t], v_cr[t], iter + 1, lr_now);
            adam_update(b_r[t], g_br[t], m_br[t], v_br[t], iter + 1, lr_now);
            adam_update(log_sd_r[t], g_lr[t], m_lr[t], v_lr[t], iter + 1, lr_now);
            log_sd_r[t] = arma::clamp(log_sd_r[t], -10.0, 5.0);
        }

        invalid_total += invalid_plus + invalid_minus + invalid_cur;
        if (invalid_cur > K_rb / 2) {
            // A self-recovery step for unstable Hamiltonian trajectories.
            theta[2 * d] = std::max(theta[2 * d] - std::log(2.0), std::log(1e-5));
            epsilon = std::exp(theta[2 * d]);
        }

        elbo_history[iter] = Lcur;
        double old_ema = elbo_ema;
        if (!std::isfinite(elbo_ema)) elbo_ema = Lcur;
        else elbo_ema = 0.95 * elbo_ema + 0.05 * Lcur;
        double dtheta = arma::norm(theta - old_theta, 2) /
            (1.0 + arma::norm(old_theta, 2));
        double delbo = std::isfinite(old_ema) ?
            std::abs(elbo_ema - old_ema) / (1.0 + std::abs(old_ema)) :
            std::numeric_limits<double>::infinity();
        criterion = std::max(dtheta, delbo);
        completed = iter + 1;

        if (Status > 0 && (iter + 1) % Status == 0) {
            Rcout << "Iteration: " << iter + 1
                  << ", auxiliary lower bound: " << elbo_ema
                  << ", epsilon: " << std::exp(theta[2 * d])
                  << ", convergence criterion: " << criterion << std::endl;
        }

        if (vb_convergence_update(criterion, Stop_Tolerance, iter,
                                  Min_Iterations, Patience, stable_count)) {
            converged = true;
            break;
        }
    }

    mu_q = theta.subvec(0, d - 1);
    log_sd_q = theta.subvec(d, 2 * d - 1);
    epsilon = std::exp(theta[2 * d]);

    // Re-estimate the final auxiliary lower bound with fresh common-random-number
    // draws.  The running EMA is useful for stopping, but it is not a clean estimate
    // of the bound at the final variational parameters.
    int K_eval = std::max(K_rb, 20);
    std::vector<arma::vec> eval_eps0(K_eval);
    std::vector<std::vector<arma::vec> > eval_momentum(
        K_eval, std::vector<arma::vec>(T_mcmc));
    for (int k = 0; k < K_eval; ++k) {
        eval_eps0[k] = arma::randn(d);
        for (int t = 0; t < T_mcmc; ++t) eval_momentum[k][t] = arma::randn(d);
    }
    int invalid_eval = 0;
    double final_bound = vb_hvi_bound(Model, Data, theta, eval_eps0, eval_momentum,
                                      a_r, c_r, b_r, log_sd_r,
                                      T_mcmc, L_leapfrog, h, invalid_eval);
    invalid_total += invalid_eval;

    arma::mat latent(n_draws, d, fill::zeros);
    arma::vec LP(n_draws, fill::zeros);
    for (int i = 0; i < n_draws; ++i) {
        arma::vec z = mu_q + arma::exp(log_sd_q) % arma::randn(d);
        for (int t = 0; t < T_mcmc; ++t) {
            arma::vec v = arma::randn(d);
            arma::vec g = vb_grad(Model, Data, z, h);
            bool valid = g.is_finite();
            arma::vec old_z = z;
            for (int l = 0; l < L_leapfrog && valid; ++l) {
                v += 0.5 * epsilon * g;
                z += epsilon * v;
                g = vb_grad(Model, Data, z, h);
                if (!g.is_finite() || !z.is_finite() || !v.is_finite()) valid = false;
                if (valid) v += 0.5 * epsilon * g;
            }
            if (!valid || !std::isfinite(vb_lp(Model, Data, z))) {
                z = old_z;
                break;
            }
        }
        vb_store_draw(Model, Data, z, i, thinned, postpred, Dev, Mon, latent, LP);
    }

    arma::vec final_mean = arma::mean(latent, 0).t();
    arma::mat final_cov;
    if (n_draws > 1) final_cov = arma::cov(latent);
    else final_cov = arma::diagmat(arma::exp(2.0 * log_sd_q));
    final_cov = vb_make_pd(final_cov, 1e-8);

    arma::mat A_r(T_mcmc, d), C_r(T_mcmc, d), B_r(T_mcmc, d), LogSD_r(T_mcmc, d);
    for (int t = 0; t < T_mcmc; ++t) {
        A_r.row(t) = a_r[t].t();
        C_r.row(t) = c_r[t].t();
        B_r.row(t) = b_r[t].t();
        LogSD_r.row(t) = log_sd_r[t].t();
    }

    arma::vec hist = elbo_history.head(std::max(completed, 1));
    return List::create(
        Named("thinned") = thinned,
        Named("postpred") = postpred,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("latent") = latent,
        Named("LP") = LP,
        Named("mean") = final_mean,
        Named("VarCov") = final_cov,
        Named("q0_mean") = mu_q,
        Named("q0_VarCov") = arma::diagmat(arma::exp(2.0 * log_sd_q)),
        Named("A_r") = A_r,
        Named("C_r") = C_r,
        Named("b_r") = B_r,
        Named("log_sd_r") = LogSD_r,
        Named("LowerBound") = final_bound,
        Named("ELBO_history") = hist,
        Named("epsilon") = epsilon,
        Named("criterion") = criterion,
        Named("n_iter") = completed,
        Named("converged") = converged,
        Named("invalid_trajectories") = invalid_total
    );
}

// [[Rcpp::export]]
List svgd(Function Model, List Data, int Iterations, int Status,
          arma::mat InitialParticles,
          int Thinning,
          arma::mat thinned, arma::mat postpred,
          arma::mat Dev, arma::mat Mon,
          double step_size = 0.01,
          double h_scale = 1.0,
          bool use_adam = true,
          double h_grad = 1e-4,
          double Stop_Tolerance = 1e-4,
          int Min_Iterations = 200,
          int Patience = 50,
          double grad_clip = 100.0) {

    RNGScope scope;
    int n = InitialParticles.n_rows;
    int d = InitialParticles.n_cols;
    arma::mat X = InitialParticles;
    arma::mat G(n, d, fill::zeros);
    arma::mat adam_m(n, d, fill::zeros);
    arma::mat adam_v(n, d, fill::zeros);
    double beta1 = 0.9, beta2 = 0.999, adam_eps = 1e-8;

    arma::vec stein_history(Iterations);
    stein_history.fill(arma::datum::nan);
    double bw = 1.0;
    double criterion = std::numeric_limits<double>::infinity();
    double criterion_ema = criterion;
    int stable_count = 0;
    int invalid_gradients = 0;
    int completed = 0;
    bool converged = false;
    double step_now = 0.0;

    for (int iter = 0; iter < Iterations; ++iter) {
        for (int i = 0; i < n; ++i) {
            arma::vec g = vb_grad(Model, Data, X.row(i).t(), h_grad);
            if (!g.is_finite()) {
                g.zeros();
                ++invalid_gradients;
            }
            double ng = arma::norm(g, 2);
            if (ng > grad_clip && ng > 0.0) g *= grad_clip / ng;
            G.row(i) = g.t();
        }
        step_now = step_size / std::sqrt(iter + 1.0);

        arma::vec x2 = arma::sum(arma::square(X), 1);
        arma::mat D = arma::repmat(x2, 1, n) +
            arma::repmat(x2.t(), n, 1) - 2.0 * X * X.t();
        D.transform([](double v) { return std::max(v, 0.0); });

        arma::vec positive = D.elem(arma::find(D > 1e-14));
        double median_sq = positive.n_elem > 0 ? arma::median(positive) : 1.0;
        bw = h_scale * median_sq / std::log((double)n + 1.0);
        bw = std::min(std::max(bw, 1e-8), 1e8);

        arma::mat K = arma::exp(-D / bw);
        arma::vec colsum = arma::sum(K, 0).t();
        arma::mat attractive = K.t() * G;
        arma::mat repulsive = (2.0 / bw) *
            (arma::diagmat(colsum) * X - K.t() * X);
        arma::mat Phi = (attractive + repulsive) / (double)n;

        arma::mat old_X = X;
        /* Adam or ordinary SVGD update occurs here */
        double update_rms = std::sqrt(
            arma::accu(arma::square(X - old_X)) /
            static_cast<double>(n * d)
        );
        double particle_scale = std::sqrt(
            arma::accu(arma::square(old_X)) /
            static_cast<double>(n * d)
        );
        criterion = update_rms / (1.0 + particle_scale);

        //criterion = std::sqrt(arma::accu(arma::square(Phi)) / (double)(n * d));
        stein_history[iter] = criterion;
        if (!std::isfinite(criterion_ema)) criterion_ema = criterion;
        else criterion_ema = 0.95 * criterion_ema + 0.05 * criterion;

        //arma::mat old_X = X;
        if (use_adam) {
            int t = iter + 1;
            adam_m = beta1 * adam_m + (1.0 - beta1) * Phi;
            adam_v = beta2 * adam_v + (1.0 - beta2) * arma::square(Phi);
            arma::mat mhat = adam_m / (1.0 - std::pow(beta1, t));
            arma::mat vhat = adam_v / (1.0 - std::pow(beta2, t));
            //X += step_size * mhat / (arma::sqrt(vhat) + adam_eps);
            X += step_now * mhat / (arma::sqrt(vhat) + adam_eps);
        } else {
            //X += step_size * Phi;
            X += step_now * Phi;
        }

        if (!X.is_finite()) {
            X = old_X;
            //step_size *= 0.5;
            step_now *= 0.5;
            adam_m.zeros();
            adam_v.zeros();
        }

        completed = iter + 1;
        if (Status > 0 && (iter + 1) % Status == 0) {
            double mean_lp = 0.0;
            for (int i = 0; i < n; ++i) mean_lp += vb_lp(Model, Data, X.row(i).t());
            mean_lp /= (double)n;
            Rcout << "Iteration: " << iter + 1
                  << ", mean LP: " << mean_lp
                  << ", Stein update RMS: " << criterion_ema
                  << ", bandwidth: " << bw << std::endl;
        }

        if (vb_convergence_update(criterion_ema, Stop_Tolerance, iter,
                                  Min_Iterations, Patience, stable_count)) {
            converged = true;
            break;
        }
    }

    arma::mat particles(n, d, fill::zeros);
    arma::mat particle_yhat(n, postpred.n_cols, fill::zeros);
    arma::mat particle_Dev(n, Dev.n_cols, fill::zeros);
    arma::mat particle_Mon(n, Mon.n_cols, fill::zeros);
    arma::vec particle_LP(n, fill::zeros);

    for (int i = 0; i < n; ++i) {
        List out = Model(wrap(X.row(i).t()), Data);
        particles.row(i) = as<arma::vec>(out["parm"]).t();
        particle_yhat.row(i) = as<arma::vec>(out["yhat"]).t();
        particle_Dev.row(i) = as<arma::vec>(out["Dev"]).t();
        particle_Mon.row(i) = as<arma::vec>(out["Monitor"]).t();
        particle_LP[i] = as<double>(out["LP"]);
    }

    // Preserve the common return fields for compatibility, but store only final
    // particles rather than the transient optimization path.
    int n_store = thinned.n_rows;
    arma::mat latent(n_store, d, fill::zeros);
    arma::vec LP(n_store, fill::zeros);
    for (int i = 0; i < n_store; ++i) {
        int j = i % n;
        thinned.row(i) = particles.row(j);
        postpred.row(i) = particle_yhat.row(j);
        Dev.row(i) = particle_Dev.row(j);
        Mon.row(i) = particle_Mon.row(j);
        latent.row(i) = X.row(j);
        LP[i] = particle_LP[j];
    }

    arma::vec hist = stein_history.head(std::max(completed, 1));
    return List::create(
        Named("thinned") = thinned,
        Named("postpred") = postpred,
        Named("Dev") = Dev,
        Named("Mon") = Mon,
        Named("latent") = latent,
        Named("LP") = LP,
        Named("particles") = particles,
        Named("latent_particles") = X,
        Named("particles_yhat") = particle_yhat,
        Named("particles_Dev") = particle_Dev,
        Named("particles_Mon") = particle_Mon,
        Named("particles_lp") = particle_LP,
        Named("bandwidth") = bw,
        Named("Stein_update_history") = hist,
        Named("criterion") = criterion_ema,
        Named("n_iter") = completed,
        Named("converged") = converged,
        Named("invalid_gradients") = invalid_gradients,
        Named("step_size") = step_now
    );
}

// ====================================THE END====================================
