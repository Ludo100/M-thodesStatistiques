InSampleStats <- function(dimension, nb_universe, returns_matrix, optimization_type, cov_matrix){
  #optimization_type: maxdiv or minvol in string
  
  number_weeks <- nrow(returns_matrix)
  
  weights <- matrix(nrow=nb_universe, ncol=dimension)
  returns <- matrix(nrow=number_weeks, ncol=nb_universe)
  
  annual_returns <- matrix(nrow=nb_universe)
  volatility <- matrix(nrow=nb_universe)
  
  for (i in 1:nb_universe){
    
    #Calculate weights of each sub-universe
    random_stocks <- sample(1:ncol(returns_matrix), size=dimension, replace=F)
    weights[i,] <- optimalPortfolio(cov_matrix[random_stocks,random_stocks], control=list(type=optimization_type, constraint = 'lo'))
    
    returns[,i] <- t(weights[i,] %*% t(returns_matrix[,random_stocks]))
    volatility[i] <- sd(returns[,i])*(sqrt(52))
    annual_returns[i] <- (cumprod(returns[,i] + 1))[number_weeks]^((1/number_weeks) * 52) - 1
  }
  
  list(returns=returns, volatility=volatility, annual_returns=annual_returns)
}


ParametricBootstrapInSample <- function(rets, risk_free, mu, sigma, n_sim, n_ptf, type, df=NA){
  #type: normal or student, if student you need to specify the degree of freedom
  #df: degree of freedom
  
  nb_stocks <- ncol(rets)
  
  #Preallocate the objects
  mu_boot <- matrix(nrow=n_sim, ncol=nb_stocks)
  sigma_boot <- array(0, c(nb_stocks, nb_stocks, n_sim))
  frontier_boot_er <- matrix(nrow=n_sim, ncol=n_ptf)
  frontier_boot_sd <- matrix(nrow=n_sim, ncol=n_ptf)
  frontier_boot_weights <- array(0, c(nb_stocks, n_ptf, n_sim))
  
  #Optimal portfolios given the simulations of the returns
  weights_minvol <- matrix(nrow=n_sim, ncol=nb_stocks)
  weights_maxdiv <- matrix(nrow=n_sim, ncol=nb_stocks)
  #Expected return for minimum volatility and maximum diversification
  mu_minvol <- rep(0, n_sim)    
  mu_maxdiv <- rep(0, n_sim)
  #Standard deviation for minimum volatility and maximum diversification
  sigma_minvol <- rep(0, n_sim)
  sigma_maxdiv <- rep(0, n_sim)
  #Sharpe Ratio for minimum volatility and maximum diversification
  sharpe_minvol <- rep(0, n_sim)
  sharpe_maxdiv <- rep(0, n_sim)
  
  set.seed(1111)
  
  #Simulate bootstrap
  if (type == "normal"){
    
    for (i in 1:n_sim) {
      
      #Simulate returns with multivariate normal
      sim_rets <- rmvnorm(mean = mu, sigma = sigma, n = nrow(rets))
      
      mu_boot[i,] <- colMeans(sim_rets)
      sigma_boot[,,i] <- cov(sim_rets)
      
      #Efficient frontiers
      frontier_boot <- efficient.frontier(mu_boot[i,], sigma_boot[,,i], nport = n_ptf, alpha.min = -0.5,
                                          alpha.max = 1.5, shorts = TRUE)
      frontier_boot_er[i,] <- frontier_boot$er
      frontier_boot_sd[i,] <- frontier_boot$sd
      frontier_boot_weights[,,i] <- t(frontier_boot$weights)
      
      weights_minvol[i,] <- optimalPortfolio(sigma_boot[,,i], control=list(type="minvol",constraint = 'lo'))
      weights_maxdiv[i,] <- optimalPortfolio(sigma_boot[,,i], control=list(type="maxdiv",constraint = 'lo'))
      returns_minvol <- weights_minvol[i,] %*% t(rets)
      returns_maxdiv <- weights_maxdiv[i,] %*% t(rets)
      
      mu_minvol[i] <- (cumprod(returns_minvol + 1))[nb_weeks]^((1/nb_weeks) * 52) - 1
      mu_maxdiv[i] <- (cumprod(returns_maxdiv + 1))[nb_weeks]^((1/nb_weeks) * 52) - 1
      sigma_minvol[i] <- sd(returns_minvol) * sqrt(52)
      sigma_maxdiv[i] <- sd(returns_maxdiv) * sqrt(52)
      sharpe_minvol[i] <- (mu_minvol[i] - risk_free) / sigma_minvol[i]
      sharpe_maxdiv[i] <- (mu_maxdiv[i] - risk_free) / sigma_maxdiv[i]
      
    }
  } else if (type == "student"){
    
    for (i in 1:n_sim) {
      
      #Simulate returns with multivariate student-t
      sim_rets <- rmvt(delta = mu, sigma = sigma, df = df, n = nrow(rets))
      
      mu_boot[i,] <- colMeans(sim_rets)
      sigma_boot[,,i] <- cov(sim_rets)
      
      #Efficient frontiers
      frontier_boot <- efficient.frontier(mu_boot[i,], sigma_boot[,,i], nport = n_ptf, alpha.min = -0.5,
                                          alpha.max = 1.5, shorts = TRUE)
      frontier_boot_er[i,] <- frontier_boot$er
      frontier_boot_sd[i,] <- frontier_boot$sd
      frontier_boot_weights[,,i] <- t(frontier_boot$weights)
      
      weights_minvol[i,] <- optimalPortfolio(sigma_boot[,,i], control=list(type="minvol",constraint = 'lo'))
      weights_maxdiv[i,] <- optimalPortfolio(sigma_boot[,,i], control=list(type="maxdiv",constraint = 'lo'))
      returns_minvol <- weights_minvol[i,] %*% t(rets)
      returns_maxdiv <- weights_maxdiv[i,] %*% t(rets)
      
      mu_minvol[i] <- (cumprod(returns_minvol + 1))[nb_weeks]^((1/nb_weeks) * 52) - 1
      mu_maxdiv[i] <- (cumprod(returns_maxdiv + 1))[nb_weeks]^((1/nb_weeks) * 52) - 1
      sigma_minvol[i] <- sd(returns_minvol) * sqrt(52)
      sigma_maxdiv[i] <- sd(returns_maxdiv) * sqrt(52)
      sharpe_minvol[i] <- (mu_minvol[i] - risk_free) / sigma_minvol[i]
      sharpe_maxdiv[i] <- (mu_maxdiv[i] - risk_free) / sigma_maxdiv[i]
    }
    
  } else{
    warning('type not supported')
  }
  
  
  mu_minvol_95 <- quantile(mu_minvol, c(0.025, 0.975))
  mu_maxdiv_95 <- quantile(mu_maxdiv, c(0.025, 0.975))
  sigma_minvol_95 <- quantile(sigma_minvol, c(0.025, 0.975))
  sigma_maxdiv_95 <- quantile(sigma_maxdiv, c(0.025, 0.975))
  sharpe_minvol_95 <- quantile(sharpe_minvol, c(0.025, 0.975))
  sharpe_maxdiv_95 <- quantile(sharpe_maxdiv, c(0.025, 0.975))
  
  
  list(frontier_boot_er=frontier_boot_er, frontier_boot_sd=frontier_boot_sd,frontier_boot_weights=frontier_boot_weights,
       weights_minvol=weights_minvol, weights_maxdiv=weights_maxdiv, returns_minvol=returns_minvol,
       returns_maxdiv=returns_maxdiv, mu_minvol=mu_minvol, mu_maxdiv=mu_maxdiv,
       sigma_minvol=sigma_minvol, sigma_maxdiv=sigma_maxdiv, sharpe_minvol=sharpe_minvol,
       sharpe_maxdiv=sharpe_maxdiv, mu_minvol_95=mu_minvol_95, mu_maxdiv_95=mu_maxdiv_95,
       sigma_minvol_95=sigma_minvol_95, sigma_maxdiv_95=sigma_maxdiv_95, sharpe_minvol_95=sharpe_minvol_95,
       sharpe_maxdiv_95=sharpe_maxdiv_95)
}

ResamplingInSample <- function(rets, risk_free, n_sim, bootstrap_id, nb_weeks, nb_stocks){
  #Function that calculates mu,sigma and sharpe of a minvol and maxdiv portfolio given the bootstrap id
  
  #Preallocate the objects
  mu_boot <- matrix(nrow=n_sim, ncol=nb_stocks)
  sigma_boot <- array(0, c(nb_stocks, nb_stocks, n_sim))
  
  #Optimal portfolios given the simulations of the returns
  weights_minvol <- matrix(nrow=n_sim, ncol=nb_stocks)
  weights_maxdiv <- matrix(nrow=n_sim, ncol=nb_stocks)
  #Expected return for minimum volatility and maximum diversification
  mu_minvol <- rep(0, n_sim)    
  mu_maxdiv <- rep(0, n_sim)
  #Standard deviation for minimum volatility and maximum diversification
  sigma_minvol <- rep(0, n_sim)
  sigma_maxdiv <- rep(0, n_sim)
  #Sharpe Ratio for minimum volatility and maximum diversification
  sharpe_minvol <- rep(0, n_sim)
  sharpe_maxdiv <- rep(0, n_sim)
  
  for (i in 1:n_sim) {
    
    mu_boot[i,] <- colMeans(as.data.frame(rets)[bootstrap_id[i,],])
    sigma_boot[,,i] <- cov(as.data.frame(rets)[bootstrap_id[i,],])
    
    weights_minvol[i,] <- optimalPortfolio(sigma_boot[,,i], control=list(type="minvol",constraint = 'lo'))
    weights_maxdiv[i,] <- optimalPortfolio(sigma_boot[,,i], control=list(type="maxdiv",constraint = 'lo'))
    returns_minvol <- weights_minvol[i,] %*% t(rets)
    returns_maxdiv <- weights_maxdiv[i,] %*% t(rets)
    
    mu_minvol[i] <- (cumprod(returns_minvol + 1))[nb_weeks]^((1/nb_weeks) * 52) - 1
    mu_maxdiv[i] <- (cumprod(returns_maxdiv + 1))[nb_weeks]^((1/nb_weeks) * 52) - 1
    sigma_minvol[i] <- sd(returns_minvol) * sqrt(52)
    sigma_maxdiv[i] <- sd(returns_maxdiv) * sqrt(52)
    sharpe_minvol[i] <- (mu_minvol[i] - risk_free) / sigma_minvol[i]
    sharpe_maxdiv[i] <- (mu_maxdiv[i] - risk_free) / sigma_maxdiv[i]
    
  }
  
  mu_minvol_95 <- quantile(mu_minvol, c(0.025, 0.975))
  mu_maxdiv_95 <- quantile(mu_maxdiv, c(0.025, 0.975))
  sigma_minvol_95 <- quantile(sigma_minvol, c(0.025, 0.975))
  sigma_maxdiv_95 <- quantile(sigma_maxdiv, c(0.025, 0.975))
  sharpe_minvol_95 <- quantile(sharpe_minvol, c(0.025, 0.975))
  sharpe_maxdiv_95 <- quantile(sharpe_maxdiv, c(0.025, 0.975))
  
  
  list(weights_minvol=weights_minvol, weights_maxdiv=weights_maxdiv, returns_minvol=returns_minvol,
       returns_maxdiv=returns_maxdiv, mu_minvol=mu_minvol, mu_maxdiv=mu_maxdiv,
       sigma_minvol=sigma_minvol, sigma_maxdiv=sigma_maxdiv, sharpe_minvol=sharpe_minvol,
       sharpe_maxdiv=sharpe_maxdiv, mu_minvol_95=mu_minvol_95, mu_maxdiv_95=mu_maxdiv_95,
       sigma_minvol_95=sigma_minvol_95, sigma_maxdiv_95=sigma_maxdiv_95, sharpe_minvol_95=sharpe_minvol_95,
       sharpe_maxdiv_95=sharpe_maxdiv_95)
}


ResamplingOutSample <- function(rets, risk_free, n_sim, bootstrap_id, nb_weeks, nb_stocks){
  #Function that calculates mu,sigma and sharpe of a minvol and maxdiv portfolio given the bootstrap id
  
  mu <- rep(0, n_sim)    
  sigma <- rep(0, n_sim)
  sharpe <- rep(0, n_sim)
  VaR5 <- rep(0, n_sim)
  
  for (i in 1:n_sim) {
    
    returns <- as.data.frame(rets)[bootstrap_id[i,],]
    
    mu[i] <- (cumprod(returns + 1))[nb_weeks]^((1/nb_weeks) * 52) - 1
    sigma[i] <- sd(returns) * sqrt(52)
    sharpe[i] <- (mu[i] - risk_free) / sigma[i]
    VaR5[i] <- quantile(returns, probs = 0.05)[1]
  }
  
  mu_95 <- quantile(mu, c(0.025, 0.975))
  sigma_95 <- quantile(sigma, c(0.025, 0.975))
  sharpe_95 <- quantile(sharpe, c(0.025, 0.975))
  VaR5_95 <- quantile(VaR5, c(0.025, 0.975))
  
  list(mu=mu, sigma=sigma, sharpe=sharpe, VaR5=VaR5,  mu_95=mu_95,  
       sigma_95=sigma_95, sharpe_95=sharpe_95,  VaR5_95=VaR5_95)
}


f_obj_maxdiv_equalw <- function(x, info) {
  # Objective function - minus the diversification ratio, as found in the class notes (p. 104)
  
  sx <- length(x) # number of assets selected
  w <- rep(1 / sx, sx) # equal-weight portfolio
  num <- crossprod(w, sqrt(diag(info$Sigma)[x]))
  den <- sqrt(tcrossprod(w, crossprod(w, info$Sigma[x, x])))
  DR <- num / den
  -DR
}

f_obj_maxdiv_optimized <- function(x, info) {
  # Objective function - minus the diversification ratio, as found in the class notes (p. 104)
  
  sx <- length(x) # number of assets selected
  w <- optimalPortfolio(info$Sigma[x,x], control=list(type="maxdiv", constraint = 'lo'))
  num <- crossprod(w, sqrt(diag(info$Sigma)[x]))
  den <- sqrt(tcrossprod(w, crossprod(w, info$Sigma[x, x])))
  DR <- num / den
  -DR
}

ParametricBootstrapOutOfSample <- function(nb_stocks, risk_free, mu, sigma, n_sim, optimization){
  #type: normal or student, if student you need to specify the degree of freedom
  #optimization: minvol or maxdiv
  #df: degree of freedom
  
  #Preallocate the objects
  #Optimal portfolios given the simulations of the returns
  weights <- matrix(nrow=n_sim, ncol=nb_stocks)
  
  set.seed(1111)
  
  #Simulate bootstrap
  
  for (i in 1:n_sim) {
    
    #Simulate returns with multivariate normal
    sim_rets <- rmvnorm(mean = mu, sigma = sigma, n = nrow(rets))
    
    sigma_boot <- cov(sim_rets)
    
    weights[i,] <- optimalPortfolio(sigma_boot, control=list(type=optimization,constraint = 'lo'))
  }
  
  weights
}


f_neighbour <- function(xc, info) {
  #Neighbour function, as found in the notes. The code was not running so we had to do some small modifications
  
  xn <- xc
  # position of assets to be changed
  # ?sample.int()
  pos <- sample.int(n = info$n_a, size = info$n_n, replace = FALSE)
  
  if (pos %in% xn){
    xn <- xn[-which(xn==pos)]
  }
  else{
    xn <- c(xn, pos)
  }
  
  sumx <- length(xn) # number of invested assets
  
  if ((sumx > info$K_sup) || (sumx < info$K_inf)) {
    out <- xc
  } else {
    out <- xn
  }
  out
}

neighbour <- function(xc, info) {
  xn <- xc
  p <- sample.int(info$n_a, size = 1L)
  xn[p] <- abs(xn[p] - 1L)
  # reject infeasible solution
  c1 <- sum(xn) >= info$nmin
  c2 <- sum(xn) <= (info$n_a - data$K_inf)
  if (c1 && c2) res <- xn else res <- xc
  as.logical(res)
}


f_x0 <- function(info){
  
  #Find random initialization
  random_size <- sample(x = info$K_inf:info$K_sup, size = 1, replace = FALSE)
  x0 <- sample(x = 1:info$n_a, size = random_size, replace = FALSE)
  x0
}
