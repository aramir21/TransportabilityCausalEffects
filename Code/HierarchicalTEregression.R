# Regression: Bayesian hierarchical models for treatment effects
rm(list = ls())
set.seed(010101)
J <- 100 # Treatments
X <- cbind(rep(1,J),matrix(rnorm(3*J),J,3))
beta <- c(1,2,-1,2); k <- length(beta)
tau <- 0.5
sig <- seq(0.1, 1.08, length.out = J) # Standard errors
TE <- X%*%beta + rnorm(J, 0, tau + sig)
plot(TE)

######## Gibbs sampler ##########
SIGMA <- diag(sig)
SIGMAi <- solve(SIGMA)

# Hyperparameters
B0 <- diag(k)
b0 <- rep(0, k)
tau0 <- 0
sigtau <- 10
  
# Gibbs functions  
PostBeta <- function(B0, b0, tau, eta){
  B0i <- solve(B0)
  Bn <- solve(B0i + t(X)%*%SIGMAi%*%X)
  bn <- Bn%*%(B0i%*%b0 + t(X)%*%SIGMAi%*%(TE - tau*eta))
  Beta <- MASS::mvrnorm(1, bn, Bn)
  return(Beta)
}

PostTau <- function(sigtau, tau0, beta, eta){
  sig2taun <- (1/sigtau^2 + as.numeric(t(eta)%*%SIGMAi%*%eta))^(-1)
  taun <- sig2taun*(tau0/sigtau^2 + as.numeric(t(eta)%*%SIGMAi%*%(TE - X%*%beta)))
  Tau <- rnorm(1, taun, sig2taun^0.5)
  return(Tau)
}

PostEta <- function(tau, beta){
  ETAn <- solve(tau^2*SIGMAi + diag(J))
  etan <- tau*ETAn%*%SIGMAi%*%(TE - X%*%beta)
  Eta <- MASS::mvrnorm(1, etan, ETAn)
  return(Eta)
}

# Posterior draws
mcmc <- 1000; burnin <- 100; tot <- mcmc+burnin; thin <- 1
PostBetas <- matrix(0, mcmc+burnin, k)
PostTaus <- rep(0, mcmc+burnin)
PostEtas <- matrix(0, mcmc+burnin, J)
Beta <- rep(0, k)
Tau <- 1
for(s in 1:tot){
  Eta <- PostEta(tau = Tau, beta = Beta) 
  Tau <- PostTau(sigtau = sigtau, tau0 = tau0, beta = Beta, eta = Eta)
  Beta <- PostBeta(B0 = B0, b0 = b0, tau = Tau, eta = Eta)
  PostEtas[s,] <- Eta
  PostTaus[s] <- Tau
  PostBetas[s,] <- Beta
}
keep <- seq((burnin+1), tot, thin)
PosteriorBetas <- PostBetas[keep,]
summary(coda::mcmc(PosteriorBetas))

PosteriorTau <- PostTaus[keep]
summary(coda::mcmc(PosteriorTau))

PosteriorEtas <- PostEtas[keep,]

# Predictive distribution
x0 <- c(1, 2, 2, 2)
S <- length(keep)
TEj <- c(sapply(1:S, function(s) {rnorm(1, mean = c(t(x0)%*%PostBetas[s, ] + PosteriorTau[s]*PosteriorEtas[s, 50]), sd = 1)}))
require(ggplot2)
dfj <- data.frame(TEj = TEj)
ggplot(dfj, aes(x = TEj)) +
  geom_histogram(binwidth = 0.2, fill = "blue", color = "black", ) +
  labs(title = "Histogram: Treatment efect x0", x = "Conditional treatment|x0", y = "Frequency") +
  geom_vline(aes(xintercept = mean(TEj)), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile(TEj, 0.025)), color = "green", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile(TEj, 0.975)), color = "green", linetype = "dashed", linewidth = 1) +
  theme_minimal()
summary(coda::mcmc(TEj))
