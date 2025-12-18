freqacc <- function(tt, aa, pp=rep(1,B), V, sw=0) {
  ## Frequentist standard deviations of Bayes estimates. For exponential
  ##    family models as in (3.1) of paper; applies both to MCMC
  ##    sampling as in section 2 or bootstrap sampling section 3.
  ## tt is vector of B simulation values of parameter of interest,
  ##    as in (2.12) of paper, either from MCMC, section 2 or
  ##    bootstrap, section 3.
  ## aa is B by p matrix of simulated row p-vectors alpha.subx(mu),(2.12),
  ##    or natural parameter vectors alpha as in exponential family
  ##    bootstrap implementation, as following (3.10).
  ## pp is B vector with elements p[i] (3.7) in bootstrap implementation
  ##    or rep(1,B) in MCMC implementation.
  ## V is mle estimate of the variance of the sufficient vector betahat
  ##    in (3.1); set to V=solve(var(aa)) if missing.
  ##
  ## Returns Bayes posterior estimate Ebayes "thetahat", (2.27) or (3.8),
  ##    its frequentist delta-method standard deviation, (2.8) or (3.11),
  ##    "cv", the internal coefficient of variation of Ebayes, (3.12)
  ##    [if using bootstrap implementation], and the usual Bayes
  ##    estimate of standard deviation of Ebayes.  B.Efron 2/8/14
  
  if (missing(V)) V <- solve(var(aa)) #only applies in bootstrap case.
  B <- length(tt)
  pp <- pp / sum(pp)
  Ebayes <- sum(pp * tt) #Bayes posterior estimate
  abar <- as.vector(pp %*% aa)
  ttt <- tt - Ebayes; aaa <- t(t(aa) - abar)
  sdbayes <- sum(pp * ttt^2)^.5
  covhat <- ttt %*% (pp * aaa) #covhat as in (3.10)
  sdfreq <- sqrt(covhat %*% V %*% t(covhat)) #as in (3.11) or Theorem 1
  
  if (var(pp) > 0) {#internal cv (3.12), from bootstrap resamples
    qq <- tt * pp; ss <- qq / mean(qq) - pp / mean(pp)
    cv <- sqrt(sum(ss^2)) / B
    v <- c(Ebayes, sdfreq, cv, sdbayes)
    names(v) <- c("Ebayes", "sdfreq", "cv", "sdbayes")
    return(v)
  }
  
  v <- c(Ebayes, sdfreq, sdbayes)
  names(v) <- c("Ebayes", "sdfreq", "sdbayes")
  v
}

