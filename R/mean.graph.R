mean.graph = function(max, distribution, parameters, showmean = FALSE){
  if (distribution == 'weibull'){
    life = seq(0, max)
    shape_para = parameters[1]
    scale_para = parameters[2]
    mx = as.numeric((scale_para/shape_para)*exp((life/scale_para)^shape_para)*Vectorize(pracma::incgam)((life/scale_para)^shape_para, 1/shape_para))
  }
  if (distribution == 'gompertz'){
    x = seq(0, max)
    a = parameters[1]
    lambda = parameters[2]
    zx = ((lambda)/a)*exp(a*x)
    mx = as.numeric(exp(zx)*(1/a)*Vectorize(pracma::incgam)(zx, 0))
  }
  if (distribution == 'gamma'){
    x = seq(0, max)
    a = parameters[1]
    lambda = parameters[2]
    sx = pgamma(x, shape = a, scale = lambda, lower.tail = FALSE)
    mx = as.numeric(((x^a)*exp(-x/lambda))/((lambda^(a-1))*gamma(a)*sx)+(lambda*a) - x)

  }
  if (distribution == 'gengamma.orig'){
    x = seq(0, max)
    upper_incomplete_gamma <- function(x,a) {
      #return (incgam(x,a))
      return (gamma(a) * pgamma(x, a, 1, lower = FALSE))
    }
    b = parameters[1]
    k = parameters[3]
    a = parameters[2]
    mx = as.numeric(a * upper_incomplete_gamma((x/a)^b, k + 1/b)/upper_incomplete_gamma((x/a)^b, k) - x)
  }
  if (distribution == 'exponential'){
    x = seq(0, max)
    lambda = parameters[1]
    mx = 1/lambda
    mx = rep(mx, max+1)
  }
  if (distribution == 'lnorm'){
    x = seq(0, max)
    mu = parameters[1]
    sigma = parameters[2]
    sx = pnorm(log(x), mean = mu, sd = sigma, lower.tail=FALSE)
    mx = as.numeric(exp(mu + (sigma^2)/2)*(pnorm(log(x), mean = (mu + sigma^2), sd = sigma, lower.tail = FALSE))/sx - x)
  }
  if (distribution == 'llogis'){
    x = seq(0, max)
    a = parameters[1]
    lambda = parameters[2]
    d = (x/lambda)^a
    sx = 1/(1+d)
    shapez = 1- (1/a)
    scalez = 1/a
    mx = as.numeric((lambda*(pi/a)/sin(pi/a))*(1+d)*(zipfR::Rbeta((1-sx), scalez, shapez, lower = FALSE)))
  }
  x = seq(0, max)
  plot(x, mx, main = 'MRL for Different Life Values', xlab = 'Life Values', ylab = 'Mean Residual Life')
  if (showmean == TRUE){
    return(mx)
  }
}



