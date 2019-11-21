InSampleStats <- function(dimension, nb_universe, returns_matrix, optimization_type){
  #optimization_type: maxdiv or minvol in string
  
  number_weeks <- nrow(returns_matrix)
  
  weights <- matrix(nrow=nb_universe, ncol=dimension)
  returns <- matrix(nrow=number_weeks, ncol=nb_universe)
  
  annual_returns <- matrix(nrow=nb_universe)
  volatility <- matrix(nrow=nb_universe)
  
  for (i in 1:nb_universe){
    
    #Calculate weights of each sub-universe
    random_stocks <- sample(1:ncol(returns_matrix), size=dimension, replace=F)
    weights[i,] <- optimalPortfolio(cov_matrix[random_stocks,random_stocks], control=list(type=optimization_type))
    
    returns[,i] <- t(weights[i,] %*% t(returns_matrix[,random_stocks]))
    volatility[i] <- sd(returns[,i])*(sqrt(52))
    annual_returns[i] <- (cumprod(returns[,i] + 1))[number_weeks]^((1/number_weeks) * 52) - 1
  }
  
  list(returns=returns, volatility=volatility, annual_returns=annual_returns)
}