YABS: Yet Another Bayesian Sampler (version 0.2.0)
=============

This package is intended to provide a flexible tool for Bayesian Modeling. It is similar to LaplacesDemon (in fact, most of our source code is inspired on this package) and fmcmc, but the MCMC algorithms are written in C++, speeding up the calculations in comparison to these packages. YABS computational efficiency does not equals JAGS' or Stan's, but it is reliable when fitting less complex models. Currently, YABS generates samples only for real-valued parameters (i.e., parameters that can have positive and negative real values).

This package should be considered experimental at this point of development. The following MCMC algorithms are implemented:

* Random-walk Metropolis (`algo = "rwm"` argument)
* Metropolis-within Gibbs (`algo = "mwg"` argument)
* Barker Proposal Metropolis (`algo = "barker"` argument)
* Oblique Hyperrectangle Slice Sampler (`algo = "ohss"` argument)

# Installation #

Using the 'remotes' package:

    install.packages("remotes")
    remotes::install_github("vthorrf/YABS")
