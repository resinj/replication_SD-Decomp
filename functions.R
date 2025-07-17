# Computation and Approximation
# Some functions that compute the distances and their decompositions

############################################################
# Functions using numerical integration via integrate():
# For the 1-Wasserstein distance and its decomposition
wd = function(qF,qG) integrate(function(x) abs(qF(x) - qG(x)), 
                               lower = 0, upper = 1, stop.on.error = FALSE)$value
wd_shift = function(qF,qG) integrate(function(x) pmax(0,pmin(qF(x/2) - qG(x/2),qF(1-x/2) - qG(1-x/2))), 
                                     lower = 0, upper = 1, stop.on.error = FALSE)$value
wd_disp = function(qF,qG) 0.5*integrate(function(x) pmax(0,qF(1-x/2) - qF(x/2) - qG(1-x/2) + qG(x/2)), 
                                        lower = 0, upper = 1, stop.on.error = FALSE)$value
wd_decomp = function(qF,qG) c(wd_shift(qF,qG),wd_shift(qG,qF),wd_disp(qF,qG),wd_disp(qG,qF))

# For the p-th power p-Wasserstein distance and its decomposition
# Defaults to the squared 2-Wasserstein distance
pwd = function(qF,qG,p = 2) integrate(function(x) (abs(qF(x) - qG(x)))^p, 
                                      lower = 0, upper = 1, stop.on.error = FALSE)$value
pwd_shift = function(qF,qG,p = 2) integrate(function(x) 
  pmax(0,pmin(sign(qF(x/2) - qG(x/2))*abs(qF(x/2) - qG(x/2))^p,
              sign(qF(1-x/2) - qG(1-x/2))*abs(qF(1-x/2) - qG(1-x/2))^p)), 
  lower = 0, upper = 1, stop.on.error = FALSE)$value
pwd_disp = function(qF,qG,p = 2) 0.5*integrate(function(x) 
  pmax(0,sign(qF(1-x/2) - qG(1-x/2))*abs(qF(1-x/2) - qG(1-x/2))^p - sign(qF(x/2) - qG(x/2))*abs(qF(x/2) - qG(x/2))^p), 
  lower = 0, upper = 1, stop.on.error = FALSE)$value
pwd_decomp = function(qF,qG,p = 2) c(pwd_shift(qF,qG,p = p),pwd_shift(qG,qF,p = p),
                                     pwd_disp(qF,qG,p = p),pwd_disp(qG,qF,p = p))

# For the Cramér distance based on the CDFs (l_2 distance)
cd_CDF = function(pF,pG,lower = -Inf,upper = Inf) integrate(function(x) (pF(x) - pG(x))^2,
                                                            lower = lower,upper = upper,stop.on.error = FALSE)$value
# For the CD decomposition
cd_shift = function(qF,qG) 0.5*integrate(function(y) 
  sapply(y, function(y) integrate(function(x,y) 
    pmax(0,pmin(qF(x/2) - qG(y/2),qF(1-x/2) - qG(1-y/2))) + pmax(0,qF(x/2) - qG(1-y/2)),
    lower = 0, upper = 1, y,stop.on.error = FALSE)$value), 
  lower = 0, upper = 1,stop.on.error = FALSE)$value
cd_disp = function(qF,qG) integrate(function(y) 
  sapply(y, function(y) integrate(function(x,y) 
    0.5*pmax(0,(qF(1-x/2) - qF(x/2)) - (qG(1-y/2) - qG(y/2))),
    lower = y, upper = 1, y,stop.on.error = FALSE)$value), 
  lower = 0, upper = 1,stop.on.error = FALSE)$value
cd_decomp = function(qF,qG) c(cd_shift(qF,qG),cd_shift(qG,qF),cd_disp(qF,qG),cd_disp(qG,qF))

################################################################################
# Functions for discrete distributions using summation 
# as detailed in Supplement S4.1
# For the 1-Wasserstein distance (area validation metric) and its decomposition
wd_discrete = function(quantiles.F, quantiles.G, # vectors of quantiles
                       levels.F,levels.G, # vectors of levels, where the quantile functions jump
                       return_decomp = TRUE){
  qF = stepfun(levels.F,quantiles.F,right = FALSE)
  qG = stepfun(levels.G,quantiles.G,right = FALSE)
  
  alphas = c(levels.F,levels.G,0.5,1)
  alphas = sort(c(alphas,1-alphas))
  
  n = length(alphas) - 1
  N = n/2
  a = alphas[-1] - alphas[1:n]
  
  wd  = sum(a*abs(qF(alphas[-(n+1)]) - qG(alphas[-(n+1)])))
  
  shift_p = 2*sum(a[1:N]*pmax(0,pmin(qF(alphas[1:N]) - qG(alphas[1:N]),qF(alphas[n:(N+1)]) - qG(alphas[n:(N+1)]))))
  shift_m = 2*sum(a[1:N]*pmax(0,pmin(qG(alphas[1:N]) - qF(alphas[1:N]),qG(alphas[n:(N+1)]) - qF(alphas[n:(N+1)]))))
  disp_p = sum(a[1:N]*pmax(0,qF(alphas[n:(N+1)]) - qG(alphas[n:(N+1)]) - qF(alphas[1:N]) + qG(alphas[1:N])))
  disp_m = sum(a[1:N]*pmax(0,qG(alphas[n:(N+1)]) - qF(alphas[n:(N+1)]) - qG(alphas[1:N]) + qF(alphas[1:N])))
  
  if(return_decomp) return(c(shift_p,shift_m,disp_p,disp_m))
  else return(wd)
}

# For the CD and its decomposition
cd_discrete = function(quantiles.F, quantiles.G, # vectors of quantiles
                       levels.F,levels.G, # vectors of levels, where the quantile functions jump
                       return_decomp = TRUE){
  qF = stepfun(levels.F,quantiles.F,right = FALSE)
  qG = stepfun(levels.G,quantiles.G,right = FALSE)
  
  alphas = c(levels.F,levels.G,0.5,1)
  alphas = sort(c(alphas,1-alphas))
  
  n = length(alphas) - 1
  N = n/2
  a = alphas[-1] - alphas[1:n]
  
  ij = expand.grid(i = 1:n,j = 1:n)
  i = ij$i
  j = ij$j
  
  cd  = sum(a[i]*a[j]*ifelse(i == j,1,2)*ifelse(i <= j, 
                                                pmax(0,qF(alphas[i]) - qG(alphas[j])),
                                                pmax(0,qG(alphas[j]) - qF(alphas[i]))))
  
  ij = expand.grid(i = 1:N,j = 1:N)
  i = ij$i
  j = ij$j
  
  shift_p = 2*sum(a[i]*a[j]*(pmax(0,pmin(qF(alphas[n+1-i]) - qG(alphas[n+1-j]),
                                         qF(alphas[i]) - qG(alphas[j]))) 
                             + pmax(0,qF(alphas[i]) - qG(alphas[n+1-j]))))
  shift_m = 2*sum(a[i]*a[j]*(pmax(0,pmin(qG(alphas[n+1-i]) - qF(alphas[n+1-j]),
                                         qG(alphas[i]) - qF(alphas[j]))) 
                             + pmax(0,qG(alphas[i]) - qF(alphas[n+1-j]))))
  disp_p = sum(a[i]*a[j]*ifelse(i == j,1,2)*ifelse(i >= j,pmax(0,qF(alphas[n+1-i]) - qG(alphas[n+1-j]) - qF(alphas[i]) + qG(alphas[j])),0))
  disp_m = sum(a[i]*a[j]*ifelse(i == j,1,2)*ifelse(i >= j,pmax(0,qG(alphas[n+1-i]) - qF(alphas[n+1-j]) - qG(alphas[i]) + qF(alphas[j])),0))
  
  
  if(return_decomp) return(c(shift_p,shift_m,disp_p,disp_m))
  else return(cd)
}

################################################################################
# Approximations
# Via linear interpolation (see Supplement S4.2)
approx_linear = function(quantiles.F, quantiles.G, # vectors of quantiles
                         betas, gammas = NULL, # quantile levels
                         support = NULL, # a vector with lower and upper bounds used for the interpolation in the tails
                         distance = "WD", # either "WD" or "CD"
                         return_decomp = TRUE){
  if(is.null(support)) support = c(min(quantiles.F,quantiles.G), max(quantiles.F,quantiles.G))
  
  qFhat = approxfun(x = c(0,betas,1), y = c(support[1], quantiles.F, support[2]))
  if(is.null(gammas)) qGhat = approxfun(x = c(0,betas,1), y = c(support[1], quantiles.G,support[2]))
  else qGhat = approxfun(x = c(0,gammas,1), y = c(support[1], quantiles.G,support[2]))
  
  if(distance == "WD"){
    if(return_decomp) return(wd_decomp(qFhat,qGhat))
    else return(wd(qFhat,qGhat))
  }
  else if(distance == "CD"){
    if(return_decomp) return(cd_decomp(qFhat,qGhat))
    else return(sum(cd_decomp(qFhat,qGhat)))
  }
  else warning("distance should be either WD or CD.")
}

# Via discretization (see Supplement S4.3)
approx_discrete = function(quantiles.F, quantiles.G, # vectors of quantiles
                              betas, gammas = NULL, # quantile levels
                              distance = "WD", # either "WD" or "CD"
                              return_decomp = TRUE){
  
  
  levels.F = (betas[-1] + betas[-length(betas)])/2
  if(is.null(gammas)) levels.G = levels.F
  else levels.G = (gammas[-1] + gammas[-length(gammas)])/2
  
  if(distance == "WD") return(wd_discrete(quantiles.F,quantiles.G,levels.F,levels.G,return_decomp))
  else if(distance == "CD") return(cd_discrete(quantiles.F,quantiles.G,levels.F,levels.G,return_decomp))
  else warning("distance should be either WD or CD.")
}

