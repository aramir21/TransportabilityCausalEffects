# Regression: Bayesian hierarchical models for treatment effects
rm(list = ls())
set.seed(010101)
J <- 100 # Treatments
mu <- 1.5
X <- cbind(rep(1,J),matrix(rnorm(3*J, mu, 1),J,3))
beta <- c(2,2.5,-0.5,2); k <- length(beta)
tau <- 0.5
sig <- 0.5
v <- 5
vj <- rgamma(J, v/2, v/2); mean(vj) # Heterogeneity
seObs <- sig/vj^0.5
plot(seObs); summary(seObs)
TE <- X%*%beta + rnorm(J, mean = 0, sd = (tau^2 + sig^2/vj)^0.5)
plot(TE)
summary(TE)

# Funtions to calculate mode and entropy of squared standard errors 
# Mode
mode_continuous <- function(x) {
  d <- density(x)  # Kernel density estimation
  mode_value <- d$x[which.max(d$y)]  # Find x where density is highest
  return(mode_value)
}
modeseObs2 <- mode_continuous(seObs^2)  # Compute mode
plot(density(seObs^2))

# Entropy
entropy <- function(x) {
  p <- table(x) / length(x)  # Compute probabilities
  ent <- -sum(p * log(p), na.rm = TRUE)  # Shannon entropy formula
  return(ent)
}
entropyseObs2 <- entropy(seObs^2)
mean(seObs^2); var(seObs^2)

HypIG <- function(hypar, modeObs, entropyObs){
  weightOpt <- entropyObs/modeObs
  shape <- hypar[1]
  rate <- hypar[2]
  if(shape <= 0 | rate <= 0){
    distEuc <- Inf
  }else{
    mode <- rate/(shape + 1)
    entropy <- shape + log(rate * gamma(shape)) - (1 + shape) * digamma(shape)
    distEuc <- (sum(weightOpt^2*(mode-modeObs)^2 + (entropy - entropyObs)^2))^0.5
  }
  return(distEuc)
}
hypar0 <- c(3, 1)
HypIG(hypar = hypar0, modeObs = modeseObs2, entropyObs = entropyseObs2)
OptimHypIG <- optim(hypar0, HypIG, method = "BFGS", control = list(maxit = 1000), modeObs = modeseObs2, entropyObs = entropyseObs2)
HypIG(hypar = OptimHypIG$par, modeObs = modeseObs2, entropyObs = entropyseObs2)
ExpIG <- rgamma(J, shape = OptimHypIG$par[1], rate = OptimHypIG$par[2])
mode_continuous(ExpIG); modeseObs2 
entropy(ExpIG); entropyseObs2

# Hyperparameters
a0 <- OptimHypIG$par[1]*2 # shape IG --> sigma2
d0 <- OptimHypIG$par[2]*2 # rate IG --> sigma2
b0 <- rep(0, k) # mean normal --> beta
B0 <- diag(k) # var normal --> beta
B0i <- solve(B0)
v <- 5 # shape = rate in heterogeneity
tau0 <- 0 # mean normal --> tau 
sigtau <- 10 # sd normal --> tau

# Posterior distributions programming the Gibbs sampling
PostBeta <- function(tau, eta, sigma2, vjs){
  SIGMA <- sigma2 * diag(1/vjs)
  SIGMAi <- solve(SIGMA)
  B0i <- solve(B0)
  Bn <- solve(B0i + t(X)%*%SIGMAi%*%X)
  bn <- Bn%*%(B0i%*%b0 + t(X)%*%SIGMAi%*%(TE - tau*eta))
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}

PostTau <- function(beta, eta, sigma2, vjs){
  SIGMA <- sigma2 * diag(1/vjs)
  SIGMAi <- solve(SIGMA)
  sig2taun <- (1/sigtau^2 + as.numeric(t(eta)%*%SIGMAi%*%eta))^(-1)
  taun <- sig2taun*(tau0/sigtau^2 + as.numeric(t(eta)%*%SIGMAi%*%(TE - X%*%beta)))
  Tau <- rnorm(1, taun, sig2taun^0.5)
  return(Tau)
}

PostEta <- function(tau, beta, sigma2, vjs){
  SIGMA <- sigma2 * diag(1/vjs)
  SIGMAi <- solve(SIGMA)
  ETAn <- solve(tau^2*SIGMAi + diag(J))
  etan <- tau*ETAn%*%SIGMAi%*%(TE - X%*%beta)
  Eta <- MASS::mvrnorm(1, etan, ETAn)
  return(Eta)
}

PostSig2 <- function(tau, beta, eta, vjs){
  an <- a0 + J
  dn <- d0 + t(TE - X%*%beta - tau*eta)%*%diag(vjs)%*%(TE - X%*%beta - tau*eta)
  Sig2 <- invgamma::rinvgamma(1, shape = an/2, rate = dn/2)
  return(Sig2)
}

PostVj <- function(sigma2, beta, tau, eta, j){
  v1n <- v + 1
  v2n <- v + sigma2^(-1)*(TE[j] - X[j,]%*%beta - tau * eta[j])^2
  taui <- rgamma(1, v1n/2, v2n/2)
  return(taui)
}

# Posterior draws
mcmc <- 10000; burnin <- 1000; tot <- mcmc+burnin; thin <- 10
PostBetas <- matrix(0, tot, k)
PostTaus <- rep(0, tot)
PostSigma2s <- rep(0, tot)
PostEtas <- matrix(0, tot, J)
PostVjs <- matrix(0, tot, J)
RegOLS <- lm(TE ~ X - 1); summary(RegOLS)
Beta <- RegOLS$coefficients
Sigma2 <- summary(RegOLS)$sigma^2/2 
Tau <- summary(RegOLS)$sigma^2/2
Eta <- rnorm(J)
pb <- winProgressBar(title = "progress bar", min = 0, max = tot, width = 300)
for(s in 1:tot){
  Vjs <- sapply(1:J, function(j){PostVj(sigma2 = Sigma2, beta = Beta, tau = Tau, eta = Eta, j)})
  Eta <- PostEta(tau = Tau, beta = Beta, sigma2 = Sigma2, vjs = Vjs) 
  Tau <- PostTau(beta = Beta, eta = Eta, sigma2 = Sigma2, vjs = Vjs)
  Sigma2 <- PostSig2(tau = Tau, beta = Beta, eta = Eta, vjs = Vjs) 
  Beta <- PostBeta(tau = Tau, eta = Eta, sigma2 = Sigma2, vjs = Vjs)
  PostVjs[s,] <- Vjs
  PostEtas[s,] <- Eta
  PostTaus[s] <- Tau
  PostSigma2s[s] <- Sigma2
  PostBetas[s,] <- Beta
  setWinProgressBar(pb, s, title=paste( round(s/tot*100, 0), "% done"))
}
close(pb)
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- coda::mcmc(PostBetas[keep,])
summary(PosteriorBetas)
plot(PosteriorBetas)
PosteriorTau2 <- (coda::mcmc(PostTaus[keep]))^2
summary(PosteriorTau2)
plot(PosteriorTau2)
PosteriorSigma <- (coda::mcmc(PostSigma2s[keep]))^0.5
summary(PosteriorSigma)
plot(PosteriorSigma)
PosteriorEtas <- PostEtas[keep,]
PosteriorVjs <- PostVjs[keep,]

# Predictive distribution
# x0 <- c(1, rep(mu, 3))
x0 <- c(1, rep(5*mu, 3))
S <- length(keep)
TEx0 <- c(sapply(1:S, function(s) {rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[s,])^0.5)}))
# Integrating out vjs
TEx0 <- rep(NA, S)
for(s in 1:S){
  TEx0s <- NULL
  for(j in 1:J){
    TEx0j <- NULL
    for(l in 1:S){
      TEx0lj <- rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
      TEx0j <- c(TEx0j, TEx0lj)
    }
    TEx0s <- c(TEx0s, mean(TEx0j))
  }
  TEx0[s] <- mean(TEx0s)
}
# Figure
require(ggplot2)
require(latex2exp) 
dfj <- data.frame(TEj = TEx0)
ggplot(dfj, aes(x = TEj)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black", ) +
  labs(title = TeX("Histogram: Treatment efect $x_0$"), x = "Conditional treatment", y = "Frequency") +
  geom_vline(aes(xintercept = mean(TEj)), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile(TEj, 0.025)), color = "green", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile(TEj, 0.975)), color = "green", linetype = "dashed", linewidth = 1) +
  theme_minimal()
summary(coda::mcmc(TEx0))

