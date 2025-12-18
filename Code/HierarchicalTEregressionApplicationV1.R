# Regression: Bayesian hierarchical models for treatment effects
rm(list = ls())
set.seed(010101)
data <- read.csv("DataSSTEJobPreferences1.csv")
data <- na.omit(data)
attach(data)
#Scale regressors
library(BoomSpikeSlab)
X <- as.matrix(scale(data[, -c(1:5)]))
niter <- 1000
y <- data$sste_jobpref

prior <- SpikeSlabPrior(cbind(1,X), y,
                        expected.model.size = ncol(X)/2, # expect nonzero predictors
                        prior.df = .01, # weaker prior than the default
                        prior.information.weight = .01,
                        diagonal.shrinkage = 0) # shrink to zero


######Estimate model########
SSBoomNew <- lm.spike(y ~ X, niter = niter, prior = prior)
#######Marginal analysis###########
Models <- SSBoomNew$beta != 0
Models <- Models[,-1]
Models[Models =="FALSE"] <- 0
Models[Models =="TRUE"] <- 1
PIP <- colMeans(SSBoomNew$beta != 0)
SummarySS <- summary(coda::mcmc(SSBoomNew$beta))
sort(PIP)
BMAglm <- BMA::bicreg(X, y, strict = FALSE, OR = 50) 
summary(BMAglm)
Reg <- lm(sste_jobpref ~ fertilityrate)
summary(Reg)

# Create the plot
ggplot(data, aes(x = fertilityrate, y = sste_jobpref)) +
  geom_point(aes(color = "Data"), size = 2) +  # Scatter plot with legend
  geom_smooth(method = "loess", aes(color = "Non-Parametric (LOESS)"), fill = "lightgreen", se = TRUE) +  # LOESS first
  geom_smooth(method = "lm", aes(color = "Linear Regression"), fill = "lightgray", se = TRUE) +  # Linear regression second
  scale_color_manual(values = c("Data" = "black", "Non-Parametric (LOESS)" = "darkgreen", "Linear Regression" = "red")) +  
  labs(
    title = "SSTE: Job Preference vs Fertility rate",
    x = "Fertility Rate",
    y = "SSTE Job Preference",
    color = ""
  ) +
  theme_minimal()
# J <- 100 # Treatments
# mu <- 1.5
# X <- cbind(rep(1,J),matrix(rnorm(3*J, mu, 1),J,3))
# beta <- c(2,2.5,-0.5,2); k <- length(beta)
# tau <- 0.5
# sig <- 0.5
# v <- 5
# vj <- rgamma(J, v/2, v/2); mean(vj) # Heterogeneity
# seObs <- sig/vj^0.5
# plot(seObs); summary(seObs)
# TE <- X%*%beta + rnorm(J, mean = 0, sd = (tau^2 + sig^2/vj)^0.5)
# plot(TE)
# summary(TE)

TE <- sste_jobpref
seObs <- sste_jobprefse
X <- cbind(1, fertilityrate)
J <- dim(X)[1]; k <- dim(X)[2]

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
summary(PosteriorTau2^0.5)
plot(PosteriorTau2^0.5)

require(ggplot2)
require(latex2exp) 
df.tau <- data.frame(tau = c(PosteriorTau2^0.5))
ggplot(df.tau, aes(x = tau)) +
  geom_density(fill = "blue", alpha = 0.3) +
  labs(title = TeX("Density Plot: $tau$"), x = TeX("$tau$"), y = "Density") +
  theme_minimal()

PosteriorSigma <- (coda::mcmc(PostSigma2s[keep]))^0.5
summary(PosteriorSigma)
plot(PosteriorSigma)
PosteriorEtas <- PostEtas[keep,]
PosteriorVjs <- PostVjs[keep,]

# Predictive distribution at Mean fertility rate
summary(fertilityrate)
x0 <- c(1, mean(fertilityrate))
S <- length(keep)
# i <- 27
# TEx0 <- c(sapply(1:S, function(s) {rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[s,i])^0.5)}))
# plot(density(TEx0))
# Integrating out vjs
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
TEx0 <- rep(NA, S)
for(s in 1:S){
  TEx0s <- NULL
  for(j in 1:J){
    TEx0j <- NULL
    for(l in 1:50){
      TEx0lj <- rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
      TEx0j <- c(TEx0j, TEx0lj)
    }
    TEx0s <- c(TEx0s, mean(TEx0j))
  }
  TEx0[s] <- mean(TEx0s)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
dfjMean <- data.frame(TEj = TEx0, x0_label = "x0 = Mean(fertilityrate)")

# Predictive distribution at Min fertility rate
x0 <- c(1, min(fertilityrate))
S <- length(keep)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
TEx0 <- rep(NA, S)
for(s in 1:S){
  TEx0s <- NULL
  for(j in 1:J){
    TEx0j <- NULL
    for(l in 1:50){
      TEx0lj <- rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
      TEx0j <- c(TEx0j, TEx0lj)
    }
    TEx0s <- c(TEx0s, mean(TEx0j))
  }
  TEx0[s] <- mean(TEx0s)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
dfjMin <- data.frame(TEj = TEx0, x0_label = "x0 = Min(fertilityrate)")

# Predictive distribution at Max fertility rate
x0 <- c(1, max(fertilityrate))
S <- length(keep)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
TEx0 <- rep(NA, S)
for(s in 1:S){
  TEx0s <- NULL
  for(j in 1:J){
    TEx0j <- NULL
    for(l in 1:50){
      TEx0lj <- rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
      TEx0j <- c(TEx0j, TEx0lj)
    }
    TEx0s <- c(TEx0s, mean(TEx0j))
  }
  TEx0[s] <- mean(TEx0s)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
dfjMax <- data.frame(TEj = TEx0, x0_label = "x0 = Max(fertilityrate)")

####### Figure #######
library(ggplot2)
library(latex2exp)
library(dplyr)
combined_df <- bind_rows(dfjMin, dfjMean, dfjMax)
# Create the plot
ggplot(combined_df, aes(x = TEj, fill = x0_label)) +
  geom_histogram(binwidth = 0.005, color = "black", alpha = 0.5, position = "identity") +
  labs(title = TeX("Histograms of treatment tffects at tifferent fertility rates"), 
       x = "Conditional treatment: Same-sex teacher job preference", 
       y = "Frequency") +
  # Add vertical lines for the mean of each group
  geom_vline(data = combined_df %>% group_by(x0_label) %>% summarise(mean_TEj = mean(TEj)),
             aes(xintercept = mean_TEj, color = x0_label), 
             linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "orange", "green")) +  # Set colors for each plot
  scale_color_manual(values = c("blue", "orange", "green")) +  # Set color for mean lines
  theme(legend.title = element_blank())  # Remove legend title

# Predictive distribution at Chad fertility rate: 6.03
x0 <- c(1, 6.6)
S <- length(keep)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
TEx0 <- rep(NA, S)
for(s in 1:S){
  TEx0s <- NULL
  for(j in 1:J){
    TEx0j <- NULL
    for(l in 1:50){
      TEx0lj <- rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
      TEx0j <- c(TEx0j, TEx0lj)
    }
    TEx0s <- c(TEx0s, mean(TEx0j))
  }
  TEx0[s] <- mean(TEx0s)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
dfjNiger <- data.frame(TEj = TEx0, x0_label = "Niger")

# Predictive distribution at Chad fertility rate: 6.03
x0 <- c(1, 0.8)
S <- length(keep)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
TEx0 <- rep(NA, S)
for(s in 1:S){
  TEx0s <- NULL
  for(j in 1:J){
    TEx0j <- NULL
    for(l in 1:50){
      TEx0lj <- rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
      TEx0j <- c(TEx0j, TEx0lj)
    }
    TEx0s <- c(TEx0s, mean(TEx0j))
  }
  TEx0[s] <- mean(TEx0s)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
dfjHK <- data.frame(TEj = TEx0, x0_label = "Hong Kong")

combined_df <- bind_rows(dfjNiger, dfjHK)
# Create the plot
ggplot(combined_df, aes(x = TEj, fill = x0_label)) +
  geom_histogram(binwidth = 0.015, color = "black", alpha = 0.5, position = "identity") +
  labs(title = TeX("Histograms of treatment effects at Hong Kong and Niger fertility rates"), 
       x = "Conditional treatment: Same-sex teacher job preference", 
       y = "Frequency") +
  # Add vertical lines for the mean of each group
  geom_vline(data = combined_df %>% group_by(x0_label) %>% summarise(mean_TEj = mean(TEj)),
             aes(xintercept = mean_TEj, color = x0_label), 
             linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = seq(-0.9, 0.3, by = 0.1)) +  # Set x-axis ticks from -0.9 to 0.3
  theme_minimal() +
  scale_fill_manual(values = c("blue", "green")) +  # Set colors for each plot
  scale_color_manual(values = c("blue", "green")) +  # Set color for mean lines
  theme(legend.title = element_blank())  # Remove legend title

# 1% students pick same job preference as their teachers
x0 <- c(1, mean(fertilityrate))
S <- length(keep)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
TEx0 <- rep(NA, S)
for(s in 1:S){
  TEx0s <- NULL
  for(j in 1:J){
    TEx0j <- NULL
    for(l in 1:50){
      TEx0lj <- rnorm(1, mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
      TEx0j <- c(TEx0j, TEx0lj)
    }
    TEx0s <- c(TEx0s, mean(TEx0j))
  }
  TEx0[s] <- mean(TEx0s)
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)
dfjMean <- data.frame(TEj = TEx0, x0_label = "x0 = Mean(fertilityrate)")
meanTE <- 0.01; sdTE <- (meanTE*(1-meanTE))^0.5
TEx0NOst <- TEx0*sdTE + meanTE
IncomeOther <- 3000 # Monthly average income other professions
IncomeTeacher <- 3500 # Monthly income teachers
BaseProp <- 0.02 # Baseline (control group) proportion students declaring want to be teachers
ControlIncome <-  IncomeOther * (1-BaseProp) + IncomeTeacher*BaseProp 
TreatedIncome <-  IncomeOther * (1-(BaseProp+TEx0NOst)) + IncomeTeacher*(BaseProp+TEx0NOst)
MoneyEffect <- TreatedIncome - ControlIncome 

dfj <- data.frame(TEj = MoneyEffect)
ggplot(dfj, aes(x = TEj)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", ) +
  labs(title = TeX("Histogram: Monetary effect conditional on mean fertility rate (1.7)"), x = "Conditional treatment (monthly USD)", y = "Frequency") +
  geom_vline(aes(xintercept = mean(TEj)), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile(TEj, 0.025)), color = "green", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile(TEj, 0.975)), color = "green", linetype = "dashed", linewidth = 1) +
  theme_minimal()
summary(MoneyEffect)
Cost <- 8 # Program cost
mean(MoneyEffect>Cost)

######### Density ############
x0 <- c(1, 1.7) # Australia
TErange <- seq(0, 0.12, 0.005) # Potential range TE
nr <- length(TErange)
TEx0Mat <- matrix(NA, S, nr)
pb <- winProgressBar(title = "progress bar", min = 0, max = S, width = 300)
for(s in 1:S){
  TEg <- NULL
  for(g in 1:nr){
    TEx0s <- NULL
    for(j in 1:J){
      TEx0j <- NULL
      for(l in 1:50){
        TEx0lj <- dnorm(TErange[g], mean = c(t(x0)%*%PostBetas[s, ]), sd = (PosteriorTau2[s] + PosteriorSigma[s]^2/PosteriorVjs[l,j])^0.5)
        TEx0j <- c(TEx0j, TEx0lj)
      }
      TEx0s <- c(TEx0s, mean(TEx0j))
    }
    TEg <- c(TEg, mean(TEx0s))
  }
  TEx0Mat[s, ] <- TEg
  setWinProgressBar(pb, s, title=paste( round(s/S*100, 0),"% done"))
}
close(pb)

library(dplyr)
library(ggplot2)
require(latex2exp)
DataDens <- tibble(t = TErange,
                   lower = apply(TEx0Mat, 2, quantile, probs = 0.025),
                   upper = apply(TEx0Mat, 2, quantile, probs = 0.975),
                   meanT = colMeans(TEx0Mat))
plot_filtering_estimates <- function(df) {
  p <- ggplot(data = df, aes(x = t)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 1, fill = "lightblue") +
    geom_line(aes(y = meanT, color = "Estimate"), linewidth = 0.5) +
    scale_color_manual(values = c("Estimate" = "blue")) +
    xlab(TeX("y")) + 
    ylab("Density") +
    labs(title = "Density: Dependent variable", color = "") + # Label for legend
    theme_minimal()
  print(p)
}
plot_filtering_estimates(DataDens)


######## Splines ###########
IdOrd <- order(fertilityrate) 
y <- sste_jobpref[IdOrd]
x <- fertilityrate[IdOrd] 
knots <- quantile(x, seq(0, 1, 0.1))
BS <- splines::bs(x, knots = knots, degree = 3, Boundary.knots = range(x), intercept = FALSE)
matplot(x, BS, type = "l", lty = 1, col = rainbow(ncol(BS)))
# Function to check if a column is constant
is_constant <- function(col) {
  return(length(unique(col)) == 1)
}
constant_columns <- apply(BS, 2, is_constant)
constant_columns
X <- BS[,-c(13:14)]
BMAglm <- BMA::bicreg(X, y, strict = FALSE, OR = 50) 
summary(BMAglm)
prior <- SpikeSlabPrior(cbind(1,X), y,
                        expected.model.size = ncol(X)/2, # expect nonzero predictors
                        prior.df = .01, # weaker prior than the default
                        prior.information.weight = .01,
                        diagonal.shrinkage = 0) # shrink to zero

######Estimate model########
SSBoomNew <- lm.spike(y ~ X, niter = niter, prior = prior)
#######Marginal analysis###########
Models <- SSBoomNew$beta != 0
Models <- Models[,-1]
Models[Models =="FALSE"] <- 0
Models[Models =="TRUE"] <- 1
PIP <- colMeans(SSBoomNew$beta != 0)
sort(PIP)
idReg <- 1# Relevant regressors, PIP > 0.5
Xnew <- cbind(1, X[, idReg])
N <- dim(Xnew)[1]
k <- dim(Xnew)[2]
# MCMC parameters
mcmc <- 5000
burnin <- 5000
tot <- mcmc + burnin
thin <- 1
# Posterior distributions using packages: MCMCpack sets the model in terms of the precision matrix
posterior  <- MCMCpack::MCMCregress(y~Xnew-1, b0=b0, B0 = B0i, c0 = a0, d0 = d0, burnin = burnin, mcmc = mcmc, thin = thin)
summary(coda::mcmc(posterior))
# Predict values with 95% credible intervals
xfit <- seq(min(x), max(x), 0.2)
H <- length(xfit)
i <- 11
idfit <- sample(1:N, H)
BSfit <- splines::bs(xfit, knots = knots, degree = 3, Boundary.knots = range(x), intercept = FALSE)
Xfit <- cbind(1, BSfit[,1]) # Relevant regressors, PIP > 0.5
Fit <- matrix(NA, mcmc, H)
# posterior[posterior > 0] <- 0
for(s in 1:mcmc){
  Fit[s,] <- Xfit%*%posterior[s,1:2]
}
# Create a data frame for ggplot
plot_data <- data.frame(x = xfit, fit = colMeans(Fit), liminf = apply(Fit, 2, quantile, 0.025), 
                        limsup = apply(Fit, 2, quantile, 0.975))

ggplot() +
  geom_line(data = plot_data, aes(x, fit), color = "blue", linewidth = 1) +  # Regression line
  geom_ribbon(data = plot_data, aes(x, ymin = liminf, ymax = limsup), fill = "blue", alpha = 0.2) +  # Confidence interval
  labs(title = "B-Spline Regression with 95% Confidence Interval",
       x = "Fertility rate",
       y = "Treatment effect: Job preference") +
  theme_minimal()
