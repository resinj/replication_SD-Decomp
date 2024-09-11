# Functions for distributions used in illustrations

################################################################################
# Normal distributions
# Std. normal
pF.std.norm = pnorm
qF.std.norm = qnorm
dF.std.norm = dnorm

# Some other normal distributions
mean.shift = 1
sd.disp = 2

# Shifted from std. normal
pF.shift.norm = function(q) pnorm(q,mean = mean.shift)
qF.shift.norm = function(p) qnorm(p,mean = mean.shift)
dF.shift.norm = function(x) dnorm(x,mean = mean.shift)

# More dispersed than std. normal
pF.disp.norm = function(q) pnorm(q,sd = sd.disp)
qF.disp.norm = function(p) qnorm(p,sd = sd.disp)
dF.disp.norm = function(x) dnorm(x,sd = sd.disp)

# Shifted and more dispersed than std. normal
pF.shift.disp.norm = function(q) pnorm(q,mean = mean.shift,sd = sd.disp)
qF.shift.disp.norm = function(p) qnorm(p,mean = mean.shift,sd = sd.disp)
dF.shift.disp.norm = function(x) dnorm(x,mean = mean.shift,sd = sd.disp)

################################################################################
# Uniform distribution
low = -1.6
up = 1.6
pF.unif = function(q) punif(q,min = low,max = up)
qF.unif = function(p) qunif(p,min = low,max = up)
dF.unif = function(x) dunif(x,min = low,max = up)

################################################################################
# Mixtures of two uniforms
pF.pwlin1 = function(q) ifelse(q < 0,0,ifelse(q < 4,0.5*q/4,ifelse(q < 6, 0.5 + (q-4)/4, 1)))
pF.pwlin2 = function(q) ifelse(q < 1,0,ifelse(q < 3,(q-1)/4,ifelse(q < 7, 0.5 + (q-3)/8, 1)))

qF.pwlin1 = function(p) ifelse(p < 0.5, 8*p, 4 + (p-0.5)*4)
qF.pwlin2 = function(p) ifelse(p < 0.5, 1 + 4*p, 3 + (p-0.5)*8)

dF.pwlin1 = function(x) ifelse(x < 0, 0, ifelse(x < 4, 0.5/4, ifelse(x < 6, 1/4, 0)))
dF.pwlin2 = function(x) ifelse(x < 1, 0, ifelse(x < 3, 1/4, ifelse(x < 7, 0.5/4, 0)))

################################################################################
# Generic functions for mixtures of uniform distributions
qF.pw.unif = function(p,bounds = c(-3,-1,0,1,3),weights = NULL){
  n = length(bounds) - 1
  if(is.null(weights)) weights = rep(1/n,n)
  if(n != length(weights)) stop("Wrong number of weights!")
  else weights = weights/sum(weights)
  cumweights = c(0,cumsum(weights))
  n_p = length(p)
  bin = rowSums(matrix(rep(cumweights,each = n_p),nrow = n_p) <= p)
  x1 = bounds[bin]
  x2 = bounds[bin+1]
  return(x1 + (x2-x1)*(p - cumweights[bin])/weights[bin])
}

# Density of a mixture of uniform distributions
dF.pw.unif = function(x,bounds = c(-3,-1,0,1,3),weights = NULL){
  n = length(bounds) - 1
  if(is.null(weights)) weights = rep(1/n,n)
  if(n != length(weights)) stop("Wrong number of weights!")
  else weights = weights/sum(weights)
  cumweights = c(0,cumsum(weights))
  n_x = length(x)
  bin = rowSums(matrix(rep(bounds,each = n_x),nrow = n_x) <= x)
  x1 = bounds[ifelse(bin == 0,1,bin)]
  x2 = bounds[bin+1]
  return(ifelse(bin == 0 | bin == n+1,0,weights[ifelse(bin == 0,1,bin)]/(x2-x1)))
}