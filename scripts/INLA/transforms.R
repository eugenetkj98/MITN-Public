# Transformed empirical logit functions, and inverses to expand gap metrics to whole number line.
# Used to improve conditioning of the regression

gap_emplogit <- function (p, epsilon = 1e-5){
  return(log((p/2 + 0.5 + epsilon)/(1 - (p/2 + 0.5) + epsilon)))
}

inv_gap_emplogit <- function(y){
  return(2*((exp(y)-0.5*exp(y)-0.5)/(1+exp(y))))
}

emplogit <- function (p, epsilon = 1e-5){
  return(log((p + epsilon)/(1 - p + epsilon)))
}

inv_emplogit <- function(y){
  return(exp(y)/(1+exp(y)))
}


## IHS version with theta adjustment to account for extreme values
# Inverse Hyperbolic sin transform
ihs <- function(x, theta){  
  return(asinh(theta * x)/theta) 
}

# Inverse of the inverse hyperbolic sin transform
inv_ihs <- function(x, theta){
    return((1/theta)*sinh(theta * x))
}

# Inverse hyperbolic sin transform-- log-likelihood
ihs_loglik <- function(theta,x){
  n <- length(x)
  xt <- ihs(x, theta)
  log.lik <- -n*log(sum((xt - mean(xt))^2))- sum(log(1+theta^2*x^2))
  return(log.lik)
}


## Alternative transforms to limit output between 0 and 1
piecewise_transform <- function(x,d){
  if ((0 < x) & (x < (1-d))){
    return(x)
  } else {
    return((1-d) + d*(1 - exp(-(x-(1-d))/d)))
  }
}

inv_piecewise_transform <- function(y,d){
  if ((0 < y) & (y < (1-d))){
    return(y)
  } else {
    return(-d * log(1 - (y-(1-d))/d) + (1-d))
  }
}

## Alternative P-transform to constrain access to [0,1]
p_transform <- function(x, mu, n = 0.5){
  if (is.na(x)){
    return(NA)
  } else {
    if (mu == 0){
      return ((x/1)^n)
    } else if (mu == 1){
      return (-((1-x)^n))
    } else {
        if (x < mu){
          return(-(((abs(x-mu))/(mu))^n))
        } else {
          return(((abs(x-mu))/((1-mu)))^n)
        }
    }
  }
}

inv_p_transform <- function(p,mu, n = 0.5){
  if (is.na(p)){
    return(NA)
  } else {
    if (p < 0){
      # return(p*mu + mu)
      return(mu - mu*(-p)^(1/n))
    } else {
      # return(p*(1-mu) + mu)
      return(mu + (1-mu)*(p^(1/n)))
    }
  }
}
