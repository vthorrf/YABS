// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

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

// ====================================THE END====================================