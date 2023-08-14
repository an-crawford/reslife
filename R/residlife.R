#' Calculating Residual Life Values
#'
#' @description Calculates residual life values over a range of values. Allows the user to specify the
#' distribution and the parameters.
#' @param values Range of values over which residual life is calculated.Usually given as a vector.
#' @param distribution Name of the distribution. Needs to be one of 'weibull', 'gompertz',
#' 'gamma', 'gengamma.orig', 'exponential', 'lnorm', or 'llogis', 'genf', or 'genf.orig'.
#' @param parameters Parameters of the survival function. Needs to be inputted
#' in order as a vector, with the name of the parameter included.
#' @param p Percentile at which to calculate residual life. Default is .5.
#' @param type Type of residual life outputted. Must be "mean", "median", "percentile",
#' or "all". Default is "mean".
#' @return The residual life for a specified sequence of values.
#' @export
#'
#' @examples residlife(values = 0:60, distribution= 'weibull', parameters = c(shape = 1.2, scale = 3))
#' residlife(values = 15:35, distribution= 'gamma', parameters =  c(shape = 1.2, rate =  1.7),
#' p = .25, type ='all')
#'
#' @references Jackson CH (2016). “flexsurv: a platform for parametric survival modeling in R.” Journal of
#' statistical software, 70.
#'
#' Poynor V (2010). “Bayesian inference for mean residual life functions in survival analysis.”
#' Masters diss., University of California, Santa Cruz.
#'
#' Prentice RL (1975). “Discrimination among some parametric models.” Biometrika, 62(3),
#' 607–614.
#'
#' Stacy EW (1962). “A generalization of the gamma distribution.” The Annals of mathematical
#' statistics, pp. 1187–1192
residlife = function(values, distribution, parameters, p = .5, type = 'mean'){
  if (distribution == 'weibull'){
    life = values
    if (names(parameters)[1]!= "shape" | names(parameters)[2]!= "scale"){
      stop("incorrect parameters entered. Parameters for weibull are shape and scale")
    }
    else{
      shape_para = parameters[1]
      scale_para = parameters[2]
      mx = as.numeric((scale_para/shape_para)*exp((life/scale_para)^shape_para)*Vectorize(pracma::incgam)((life/scale_para)^shape_para, 1/shape_para))
      S_t = exp(-((life/scale_para)^shape_para))
      px = function(p){
        pc = (1-p)*S_t
        px = as.numeric(qweibull(pc, shape = shape_para, scale = scale_para, lower.tail = FALSE) - life)
      }
    }
  }
  if (distribution == 'gompertz'){
    x = values
    if (names(parameters)[1]!= "shape" | !(names(parameters)[2] %in% c("scale", "rate"))){
      stop("incorrect parameters entered. Parameters for gompertz are shape and scale/rate")
    }
    else{
      a = parameters[1]
      lambda = parameters[2]
      zx = ((lambda)/a)*exp(a*x)
      sx = exp(((lambda)/a)*(1-exp(a*x)))
      mx = as.numeric(exp(zx)*(1/a)*Vectorize(pracma::incgam)(zx, 0))
      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qgompertz(pc, shape = a, rate = lambda, lower.tail = FALSE) - x)
      }
    }
  }
  if (distribution == 'gamma'){
    x = values
    if (names(parameters)[1]!= "shape" | !(names(parameters)[2] %in% c("scale", "rate"))){
      stop("incorrect parameters entered. Parameters for gamma are shape and scale/rate")
    }
    else{
      if (names(parameters)[2] == 'rate'){
        lambda = 1/parameters[2]
      }
      else{
        lambda = parameters[2]
      }
      a = parameters[1]
      sx = pgamma(x, shape = a, scale = lambda, lower.tail = FALSE)
      mx = as.numeric(((x^a)*exp(-x/lambda))/((lambda^(a-1))*gamma(a)*sx)+(lambda*a) - x)
      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qgamma(pc, shape = a, scale = lambda, lower.tail = FALSE) - x)
      }
    }
  }
  if (distribution == 'gengamma.orig'){
    x = values
    upper_incomplete_gamma <- function(x,a) {
      #return (incgam(x,a))
      return (gamma(a) * pgamma(x, a, 1, lower.tail = FALSE))
    }
    if (names(parameters)[1]!= "shape" | names(parameters)[2]!= "scale" | names(parameters)[3]!= 'k'){
      stop("incorrect parameters entered. Parameters for gengamma.orig are shape, scale, and k")
    }
    else{
      b = parameters[1]
      k = parameters[3]
      a = parameters[2]
      mx = as.numeric(a * upper_incomplete_gamma((x/a)^b, k + 1/b)/upper_incomplete_gamma((x/a)^b, k) - x)
      sx = pgengamma.orig(x, shape = b, scale = a, k = k, lower.tail = FALSE)
      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qgengamma.orig(pc, shape = b, scale = a, k = k, lower.tail = FALSE) - x)
      }
    }
  }
  if (distribution == 'gengamma'){
    x = values
    upper_incomplete_gamma <- function(x,a) {
      return (gamma(a) * pgamma(x, a, 1, lower.tail = FALSE))
    }
    if (names(parameters)[1]!= "mu" | names(parameters)[2]!= "sigma" | names(parameters)[3]!= 'Q'){
      stop("incorrect parameters entered. Parameters for gengamma are mu, sigma, and Q")
    }
    else{
      mu = parameters[1]
      sigma = parameters[2]
      Q = parameters[3]
      k = 1 / (Q^2)
      b = abs(Q)/(sigma)
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))
      sx = pgengamma(x,mu = mu, sigma = sigma, Q = Q, lower.tail = FALSE)
      mx = as.numeric(a * upper_incomplete_gamma((x/a)^b, k + 1/b)/upper_incomplete_gamma((x/a)^b, k) - x)

      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qgengamma(pc,mu = mu, sigma = sigma,Q = Q, lower.tail = FALSE) - x)
      }

    }
  }
  if (distribution == 'exponential'){
    x = values
    if (names(parameters)[1]!= 'rate'){
      stop("incorrect parameter entered. The parameter for exponential is rate")
    }
    else{
      lambda = parameters[1]
      mx = 1/lambda
      mx = rep(mx, (max-min))
      sx = exp(-(1/lambda)*x)
      px = function(p){
        pc = (1-p)*sx
        px = qexp(pc, 1/lambda, lower.tail = FALSE)
      }
    }
  }
  if (distribution == 'lnorm'){
    x = values
    if (names(parameters)[1]!= "meanlog" | names(parameters)[2]!= "sdlog"){
      stop("incorrect parameters entered. The parameters for lnorm are meanlog and sdlog")
    }
    else{
      mu = parameters[1]
      sigma = parameters[2]
      sx = pnorm(log(x), mean = mu, sd = sigma, lower.tail=FALSE)
      mx = as.numeric(exp(mu + (sigma^2)/2)*(pnorm(log(x), mean = (mu + sigma^2), sd = sigma, lower.tail = FALSE))/sx - x)
      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qlnorm(pc, meanlog = mu, sdlog = sigma, lower.tail = FALSE)-x)
      }
    }
  }
  if (distribution == 'llogis'){
    regbeta = function (x, a, b, lower = TRUE, log = !missing(base), base = exp(1))
    {
      if (log) {
        if (missing(base)) {
          pbeta(x, shape1 = a, shape2 = b, lower.tail = lower,
                log.p = TRUE)
        }
        else {
          pbeta(x, shape1 = a, shape2 = b, lower.tail = lower,
                log.p = TRUE)/log(base)
        }
      }
      else {
        pbeta(x, shape1 = a, shape2 = b, lower.tail = lower,
              log.p = FALSE)
      }
    }
    x = values
    if (names(parameters)[1]!= "shape" | names(parameters)[2]!= "scale"){
      stop("incorrect parameters entered. The parameters fir llogis are shape and scale.")
    }
    else{
      a = parameters[1]
      lambda = parameters[2]
      d = (x/lambda)^a
      sx = 1/(1+d)
      shapez = 1- (1/a)
      scalez = 1/a
      mx = as.numeric((lambda*(pi/a)/sin(pi/a))*(1+d)*(regbeta((1-sx), scalez, shapez, lower = FALSE)))
      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qllogis(pc, shape = a, scale = (lambda), lower.tail = FALSE) - x)
      }
    }
  }
  if (distribution == 'genf'){
    x = values
    if(names(parameters)[1]!= 'mu' | names(parameters)[2]!= 'sigma' | names(parameters)[3]!= 'Q' | names(parameters)[4]!= 'P'){
      stop("incorrect parameters entered. The parameters for genf are mu, sigma, Q, and P")
    }
    else{
      mu = parameters[1]
      sigma = parameters[2]
      Q = parameters[3]
      P = parameters[4]
      delta = (Q^2 + 2*P)^(1/2)
      sigma = sigma/delta
      m1 = 2*(Q^2 + 2*P + Q*delta)^(-1)
      m2 = 2*(Q^2 + 2*P - Q*delta)^(-1)
      sx <- pgenf.orig(x, mu =  mu, sigma = sigma, s1 = m1, s2 =  m2, lower.tail = FALSE)
      C = exp(-mu/sigma)*(m1/m2)*(x^(1/sigma))
      part1 = (exp(mu)*((m2/m1)^sigma)*C^(-m2+sigma))/(beta(m1,m2))
      part2 = gamma(m2-sigma)/gamma(m2-sigma+1)
      part3 = Gauss2F1b(m1+m2, m2-sigma, m2-sigma+1, -1/C)
      num_integral = part1 * part2 * part3
      f0 <- function(x,mu) {
        return(x*dgenf.orig(x,mu = mu, sigma = sigma, s1 = m1,  s2 = m2))
      }
      integral_fun <- function(mu) integrate(f0, 0, Inf,mu)$value
      vec_integral_fun <- Vectorize(integral_fun)
      mx <- ifelse(C == 0, vec_integral_fun(mu), num_integral/sx - x)
      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qgenf.orig(pc,   mu =  mu, sigma = sigma,   s1 = m1, s2 =  m2,lower.tail = FALSE) - x)
      }
    }
  }
  if (distribution == 'genf.orig'){
    x = values
    if(names(parameters)[1]!= 'mu' | names(parameters)[2]!= 'sigma' | names(parameters)[3]!= 'm1' | names(parameters)[4]!= 'm2'){
      stop("incorrect parameters entered. The parameters for genf.orig are mu, sigma, m1, and m2.")
    }
    else{
      sx = pgenf.orig(x, mu = mu,  sigma = sigma, s1 = m1, s2 = m2, lower.tail = FALSE)
      C = exp(-mu/sigma)*(m1/m2)*(x^(1/sigma))
      part1 = (exp(mu)*((m2/m1)^sigma)*C^(-m2+sigma))/(beta(m1,m2))
      part2 = gamma(m2-sigma)/gamma(m2-sigma+1)
      part3 = Gauss2F1b(m1+m2, m2-sigma, m2-sigma+1, -1/C)
      num_integral = part1 * part2 * part3
      f0 <- function(x,mu) {
        return(x*dgenf.orig(x,mu = mu,sigma = sigma, s1 = m1,s2 = m2))
      }
      integral_fun <- function(mu) integrate(f0, 0, Inf,mu)$value
      vec_integral_fun <- Vectorize(integral_fun)
      mx <- ifelse(C == 0, vec_integral_fun(mu), num_integral/sx - x)
      px = function(p){
        pc = (1-p)*sx
        px = as.numeric(qgenf.orig(pc,  mu = mu, sigma = sigma, s1 = m1, s2 = m2,lower.tail = FALSE) - x)
      }
    }
  }
  if (type=='mean'){
    return(c(mx))
  }
  if (type== 'median'){
    return(c(px(.5)))
  }
  if (type == 'percentile'){
    return(c(px(p)))
  }
  if (type == 'all'){
    return(data.frame(mean = mx, median = px(.5), percentile = px(p)))
  }
  else{
    return('invalid type')
  }

}




