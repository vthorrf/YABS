###===--- Test File ---===###

### Start
rm(list=ls())
dev.off()
cat("\014")

### Packages and functions====
require(YABS)
require(compiler)

### Random data====
seed <- 1234; N <- 200; V <- 2
set.seed(seed)
X <- MASS::mvrnorm(N, rep(0,V), diag(V))
betas <- runif(V+1, .3, .7)
y <- rnorm(N, cbind(1,X) %*% betas, 1)

### DATA LIST====
parm.names <- c(paste("beta",0:V,sep=""),"sigma")
pos.beta   <- grep("beta", parm.names)
pos.sigma  <- grep("sigma", parm.names)
PGF <- function(Data) {
  beta  <- rnorm(ncol(Data$X)+1)
  sigma <- rnorm(1)
  return(c(beta, sigma))
}
Data <- list( X=X, y=y, n=N, parm.names=parm.names, pos.beta=pos.beta,
              pos.sigma=pos.sigma, PGF=PGF, mon.names="LP" )
Initial.Values <- PGF(Data)

### MODEL====
Model <- function(parm, Data){
  ### Parameters
  beta  <- parm[Data$pos.beta]
  sigma <- exp(parm[Data$pos.sigma])
  
  ### Log-Prior
  beta.prior  <- sum(dnorm(beta, 0, 1, log=T))
  sigma.prior <- sum(dgamma(sigma, 1e-2, 1e-2, log=T))
  Lp <- beta.prior + sigma.prior
  
  ### Log-Likelihood
  mu <- cbind(1,Data$X) %*% beta
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  
  ### Log-Posterior
  LP <- LL + Lp
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP, yhat=rnorm(Data$n, mu, sigma), parm=parm)
  return(Modelout)
}
Model <- compiler::cmpfun(Model)

### FIT====
fit1 <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
             burnin=100, status=110, thinning=1, algo="harmwg")
fit2 <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
             burnin=100, status=110, thinning=1, algo="harm")
fit3 <- MCMC(Model, Data, Initial.Values=NULL, iterations=1000,
             burnin=100, status=110, thinning=1, algo="gcharm")
plot.ts(fit1$posterior)
plot.ts(fit2$posterior)
plot.ts(fit3$posterior)
fit1
fit2
fit3

### JAGS====
require(jagsUI)
p <- V + 1 # Quantidade de par?metros beta
dataList = list( y=y, X=cbind(1,X), N=N, p=p )

## Definir modelo
RLST <- " ### Regressão linear simples tradicional
model {
  # Priors
  for (i in 1:p) {
    betas[i] ~ dnorm(0,1)
  }
  sigma      ~ dgamma(1e-2,1e-2)
  
  # Verossimilhan?a
  for(i in 1:N) {
    for(k in 1:p) {
      mu[i,k] <- betas[k]*X[i,k]
    }
    y[i] ~ dnorm(sum(mu[i,1:p]),sigma)
  }
}"

## Run the chains
# Name the parameters to be monitored
params <- c("betas","sigma")
# Random initial values
inits <- function(){ 
  sigma <- rgamma(1,1e-2,1e-2)
  betas <- rnorm(p)
  return(list(sigma=sigma, betas=betas))
}
# Define some MCMC parameters for JAGS
nthin    = 1    # How Much Thinning?
nchains  = 1    # How Many Chains?
nburnin  = 100  # How Many Burn-in Samples?
nsamples = 1100 # How Many Recorded Samples?
nadapt   = 2000 # How Many adaptation Samples?
set.seed(666)
model = textConnection(RLST)
UIsamples <- jagsUI::jags(dataList, inits=inits, params, model.file=model,
                          n.chains=nchains, n.adapt=nadapt, n.iter=nsamples,
                          n.burnin=nburnin, n.thin=nthin, DIC=T)

## Compare run time
fit1$mcmc.info$elapsed.mins
fit2$mcmc.info$elapsed.mins
fit3$mcmc.info$elapsed.mins
UIsamples$mcmc.info$elapsed.mins

