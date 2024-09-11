# Some functions to be used in our studies

############################################################
# Computation of the 1-Wasserstein distance and its decomposition using integrate()
wd = function(qF,qG) integrate(function(x) abs(qF(x) - qG(x)), 
                               lower = 0, upper = 1, stop.on.error = FALSE)$value
wd_shift = function(qF,qG) integrate(function(x) pmax(0,pmin(qF(x/2) - qG(x/2),qF(1-x/2) - qG(1-x/2))), 
                                     lower = 0, upper = 1, stop.on.error = FALSE)$value
wd_disp = function(qF,qG) 0.5*integrate(function(x) pmax(0,qF(1-x/2) - qF(x/2) - qG(1-x/2) + qG(x/2)), 
                                        lower = 0, upper = 1, stop.on.error = FALSE)$value
wd_decomp = function(qF,qG) c(wd_shift(qF,qG),wd_shift(qG,qF),wd_disp(qF,qG),wd_disp(qG,qF))

# Computation of the p-th power p-Wasserstein distance and its decomposition using integrate()
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

# Computation of the Cramér distance based on the CDFs (l_2 distance) using integrate()
cd_CDF = function(pF,pG,lower = -Inf,upper = Inf) integrate(function(x) (pF(x) - pG(x))^2,
                                                            lower = lower,upper = upper,stop.on.error = FALSE)$value
# Computation of the CD decomposition using integrate()
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
# An approximation of the Cramér distance and its decomposition
# based on a finite number of given quantiles
# as described in Supplement S4.1
cd_approx = function(quantiles.F, quantiles.G, # vectors of quantiles
                     alphas, betas = NULL, # quantile levels
                     return_decomp = TRUE){
  K = length(alphas)
  # Check required symmetries: Quantile levels (alphas) should be symmetric around 0.5 
  # for the decomposition to make sense.
  if(any(round(alphas + alphas[K:1],10) != 1))
    warning("Quantile levels do not bound central prediction intervals.")
  
  if(length(quantiles.F) != length(alphas)) 
    stop("Number of quantiles of F (quantiles.F) does not match number of quantile levels
         (alphas).")
  if(is.null(betas)){
    betas = alphas
    print("Quantile levels of G (betas) set to match quantile levels of F (alphas).")
  }
  L = length(betas)
  if(any(round(betas + betas[L:1],10) != 1)) warning("Quantile levels (betas) not symmetric.")
  if(length(quantiles.G) != length(betas)) 
    stop("Number of quantiles of G (quantiles.G) does not match number of quantile levels
         (betas/alphas).")
  
  alphas_ext = sort(unique(round(c(0,alphas,1-alphas,1),digits = 10)))
  betas_ext = sort(unique(round(c(0,betas,1-betas,1),digits = 10)))
  
  if(return_decomp){
    # Approximate components
    integrand_comps = function(i,j){
      weight_a = (alphas_ext[i+2] - alphas_ext[i])/2
      weight_b = (betas_ext[j+2] - betas_ext[j])/2
      
      lF = quantiles.F[i]
      uF = quantiles.F[K+1 - i]
      lG = quantiles.G[j]
      uG = quantiles.G[L+1 - j]
      
      return(ifelse(round(alphas[i],10) == 0.5 || round(betas[j],10) == 0.5, 1, 2)*
               # correction factor for the median times factor 2
               weight_a*weight_b*
               c(SFG = pmax(0,pmin(lF - lG, uF - uG) + pmax(0,lF - uG)),
                 SGF = pmax(0,pmin(lG - lF, uG - uF) + pmax(0,lG - uF)),
                 ifelse(alphas[i] == betas[j],0.5,1)* # correction at equal levels
                   c(DFG = ifelse(betas[j] <= alphas[i], pmax(0, (uF - lF) - (uG - lG)), 0),
                     DGF = ifelse(alphas[i] <= betas[j], pmax(0, (uG - lG) - (uF - lF)), 0))))
    }
    
    comps = rowSums(apply(expand.grid(1:ceiling(K/2),1:ceiling(L/2)), 1,
                          function(x) integrand_comps(x[1],x[2])))
    return(comps)
  }
  else{
    # Approximate CD
    integrand_comps = function(i,j){
      weight_a = (alphas_ext[i+2] - alphas_ext[i])/2
      weight_b = (betas_ext[j+2] - betas_ext[j])/2
      
      lF = quantiles.F[i]
      uF = quantiles.F[K+1 - i]
      lG = quantiles.G[j]
      uG = quantiles.G[L+1 - j]
      
      return(ifelse(round(alphas[i],10) == 0.5 || round(betas[j],10) == 0.5, 1, 2)*
               # correction factor for the median times factor 2
               weight_a*weight_b*
               (ifelse(alphas[i] == betas[j],0.5,1)* # correction at equal levels
                  (ifelse(sign(alphas[i] - betas[j]) != sign(lF - lG),abs(lF - lG),0) +
                     ifelse(sign(alphas[K+1 - i] - betas[L+1 - j]) != sign(uF - uG),abs(uF - uG),0)) +
                  ifelse(sign(alphas[i] - betas[L+1 - j]) != sign(lF - uG),abs(lF - uG),0) +
                  ifelse(sign(alphas[K+1 - i] - betas[j]) != sign(uF - lG),abs(uF - lG),0)))
    }
    
    cd = sum(apply(expand.grid(1:ceiling(K/2),1:ceiling(L/2)), 1,
                   function(x) integrand_comps(x[1],x[2])))
    return(cd)
  }
}









