YABS: Yet Another Bayesian Sampler (version 0.1.0)
=============

This package is intended to provide a flexible tool for Bayesian Modeling. It is similar to LaplacesDemon and fmcmc, but the MCMC algorithms are written in C++, speeding up the calculations in comparison to these packages.

This package should be considered experimental in this point of development. The following MCMC algorithms are implemented::

* Hit-and-Run Metropolis (`algo = "harm"` argument)
* Hit-and-Run Metropolis-within Gibbs (`algo = "harmwg"` argument)

# Installation #

Using the 'remotes' package:

    install.packages("remotes")
    remotes::install_github("vthorrf/YABS")
