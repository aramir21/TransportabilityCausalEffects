# Bayesian hierarchical models for treatment effects

TE <- c(28, 8, -3, 7, -1, 1, 18, 12)
mean(TE)
StEr <- c(15, 10, 16, 11, 9, 11, 10, 18)
J <-length(TE)

# Posterior tau
PostTau <- function(tau){
  Vmui <- sum((StEr^2+tau^2)^-1)
  muhat <- sum((1/(StEr^2+tau^2))*TE)/sum((1/(StEr^2+tau^2)))
  term <- NULL
  for(j in 1:J){
    termj <- (StEr[j]^2+tau^2)^(-1/2)*exp(-(TE[j]-muhat)^2/(2*(StEr[j]^2+tau^2)))
    term <- c(term, termj)
  }
  denstau <- Vmui^(-0.5)*prod(term)
  return(denstau)
}
tau0 <- seq(0.1,30,0.1)
ProbTau <- sapply(tau0, PostTau)
plot(tau0, ProbTau, type = "l")

PostMu <- function(tau){
  Vmui <- sum((StEr^2+tau^2)^-1)
  muhat <- sum((1/(StEr^2+tau^2))*TE)/sum((1/(StEr^2+tau^2)))
  mu <- rnorm(1, muhat, (1/Vmui)^0.5)
  return(mu)
}
tau <- 5
mu <- PostMu(tau)

PostThetaj <- function(mu, tau, j){
  thetahatj <- (TE[j]/StEr[j]^2+mu/tau^2)/(1/StEr[j]^2+1/tau^2)
  Vj <- 1/(1/StEr[j]^2+1/tau^2)
  thetaj <- rnorm(1, thetahatj, Vj^0.5)
  return(thetaj)
}
PostThetaj(mu, tau, 7)

# Gibbs sampler
S <- 10000
THETA <- matrix(NA, S, J)
tau0 <- seq(0.05,30,0.5)
ProbTau <- sapply(tau0, PostTau)
tau <- sample(tau0, S, prob = ProbTau, replace = TRUE)
mu <- sapply(tau, PostMu)
for(j in 1:J){
  thetaj <- sapply(1:S, function(s){PostThetaj(mu = mu[s], tau = tau[s], j = j)})
  THETA[,j] <- thetaj
}

summary(coda::mcmc(tau))
summary(coda::mcmc(mu))
summary(coda::mcmc(THETA))

# Posterior distribution of thetaj given tau, and averaging over mu
THETAtau <- matrix(NA, length(tau0), J)
for (s in 1:length(tau0)){
  mus <- replicate(20000, PostMu(tau0[s]))
  for(j in 1:J){
    thetajtau <- mean(sapply(1:S, function(l){PostThetaj(mu = mus[l], tau = tau0[s], j = j)}), na.rm = TRUE)
    THETAtau[s,j] <- thetajtau
  }
}

plot(tau0, THETAtau[,1], type = "l")

# Predictive distribution
j <- 3
TEj <- c(sapply(1:S, function(s) {rnorm(1, mean = THETA[s,j], sd = StEr[j])}))
require(ggplot2)
dfj <- data.frame(TEj = TEj)
ggplot(dfj, aes(x = TEj)) +
  geom_histogram(binwidth = 2, fill = "blue", color = "black", ) +
  labs(title = "Histogram: Treatment efect j", x = "Values", y = "Frequency") +
  theme_minimal()
summary(TEj)
