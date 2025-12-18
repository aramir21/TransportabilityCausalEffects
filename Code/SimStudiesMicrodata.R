set.seed(010101)
J <- 5 # Studies
N <- c(500, 256, 145, 760, 678) # Sample size microdata
k <- 3
j <- unlist(sapply(1:J, function(l){rep(l, N[l])}))
i <- unlist(sapply(1:J, function(l){1:N[l]}))
d1 <- j == 1; d2 <- j == 2; d3 <- j == 3; d4 <- j == 4; d5 <- j == 5
NJ <- length(i)
X <- matrix(rnorm(NJ*k), NJ, k)
w1 <- unlist(sapply(1:J, function(l){rep(rnorm(1), N[l])})) 
w2 <- unlist(sapply(1:J, function(l){rep(rnorm(1), N[l])}))
w3 <- unlist(sapply(1:J, function(l){rep(rnorm(1), N[l])})) 
w4 <- unlist(sapply(1:J, function(l){rep(rnorm(1), N[l])}))
w5 <- unlist(sapply(1:J, function(l){rep(rnorm(1), N[l])})) 
w6 <- unlist(sapply(1:J, function(l){rep(rnorm(1), N[l])}))
W <- cbind(w1, w2)
Tr <- rbinom(NJ, 1, 0.3)
XWT <- cbind(1, X, w1, w2, Tr)
B <- rep(1, k + 4)
e <- rnorm(NJ) 
y <- XWT%*%B + e
# We can introduce within invariable regressors as studies we have. 
# Within invariable regressors play exactly the same role as fixed studies effects, but they can give the right signal
Reg <- lm(y ~ X + Tr + w1 + w2 + w3 + w4 + w5 + w6)
summary(Reg)
Data <- data.frame(id = j, y = y, x1 = X[,1], x2 = X[,2], w1 = w1, w2 = w2, w3 = w3, w4 = w4, w5 = w5, w6 = w6, Tr = Tr)
K1 <- 6; K2 <- 1
b0 <- rep(1, K1); B0 <- diag(K1)
r0 <- 5; R0 <- diag(K2)
a0 <- 0.001; d0 <- 0.001
Resultshreg <- MCMCpack::MCMChregress(fixed = y ~ x1 + x2 + Tr + w1 + w2, random = ~1, 
                                      group = "id", data = Data, burnin = 5000, mcmc = 10000, thin = 1,
                                      r = r0, R = R0, nu = a0, delta = d0)
Betas <- Resultshreg[["mcmc"]][,1:K1]
# Sigma2RanEff <- Resultshreg[["mcmc"]][,54]
# Sigma2 <- Resultshreg[["mcmc"]][,55]
summary(Betas)
# summary(Sigma2RanEff)
# summary(Sigma2)
# summary(Sigma2RanEff/(Sigma2RanEff+Sigma2))
posterior  <- MCMCpack::MCMCregress(y ~ x1 + x2 + Tr + w1 + w2, data = Data)
summary(posterior)
