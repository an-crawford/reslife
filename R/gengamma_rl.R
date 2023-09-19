####################
##### GenGamma #####
####################

############################
#### Author: Zekai Wang ####
############################

#the parameters relationship between gengamma and gengamma.orig is:
#dgengamma.orig(x, shape=abs(Q)/(sigma), scale=exp(mu - (log(Q^(-2))*sigma)/(abs(Q))), k=1/(Q^2)) == dgengamma(x, mu=mu, sigma=sigma, Q=abs(Q))
# ------------------------------------------
#dgengamma.orig(x, shape=b, scale=a, k) == dgengamma(x, mu=log(a) + log(k)/b, sigma=1/(b*sqrt(k)), Q=1/sqrt(k))
######### example #######
#library('flexsurv')
#fs1 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc ,dist = "genGamma")
#fs2 <- flexsurvreg(Surv(recyrs, censrec) ~ 1, data = bc ,dist = "genGamma.orig")

#mu <- as.numeric(fs1$coefficients[1])
#sigma <- exp(as.numeric(fs1$coefficients[2]))
#Q <- as.numeric(fs1$coefficients[3])
#x = 4
#dgengamma.orig(x, shape=abs(Q)/(sigma), scale=exp(mu - (log(Q^(-2))*sigma)/(abs(Q))), k=1/(Q^2))
#dgengamma(x, mu=mu, sigma=sigma, Q=abs(Q))

#b <- exp(as.numeric(fs2$coefficients[1]))
#a <- exp(as.numeric(fs2$coefficients[2]))
#k <- exp(as.numeric(fs2$coefficients[3]))

#dgengamma.orig(x, shape=b, scale=a, k)
#dgengamma(x, mu=log(a) + log(k)/b, sigma=1/(b*sqrt(k)), Q=1/sqrt(k))
# ---------------------------------------------
# ---------------------------------------------
# Accuracy Test #
#fs1 <- flexsurvreg(Surv(years, status) ~ 1, data = bosms3 ,dist = "genGamma")
#m0 <- gengamma(fs1, 0, p=.5)

#f0 <- function(x) {
#  return(x*dgengamma(x,mu = fs1$coefficients[1],sigma = exp(fs1$coefficients[2]),Q = fs1$coefficients[3]))
#}
#mean <- integrate(f0, 0, Inf)$value
#m0 == mean#
# ---------------------------------------------
# ---------------------------------------------
#fs1 <- flexsurvreg(Surv(years, status) ~ 1, data = bosms3 ,dist = "genGamma")
#mu <- as.numeric(fs1$coefficients[1])
#sigma <- exp(as.numeric(fs1$coefficients[2]))
#Q <- as.numeric(fs1$coefficients[3])
##f1 <- function(x) {
#return(pgengamma(x,mu = mu,sigma = sigma,Q = Q, lower.tail = FALSE))
#}
#mrl <- gengamma(fs1, 4, p=.5, type ='mean')
#ni <- integrate(f1, 4, Inf)$value/pgengamma(4, mu = mu,sigma = sigma,Q = Q, lower.tail = FALSE)
#ni == mrl


upper_incomplete_gamma <- function(x,a) {

  #return (incgam(x,a))

  return (gamma(a) * pgamma(x, a, 1, lower.tail = FALSE))

}

gengamma_rl = function(fsroutput, x, p=.5, type = 'all', newdata = data.frame()){

  if (length(newdata)!=0){

    if (length(newdata) == 1){

      if(fsroutput$covdata$covnames != colnames(newdata)){
        stop('Wrong columns in inputted data')
      }

    }

    else{

      names = fsroutput$covdata$covnames

      newdata= newdata[,c(names)]

    }

  }

  sigma = exp(as.numeric(fsroutput$coefficients[2]))
  Q = as.numeric(fsroutput$coefficients[3])

  k = 1 / (Q^2)
  b = abs(Q)/(sigma)

  if (length(newdata) == 0){

    if (fsroutput$ncovs == 0) {

      mu = as.numeric(fsroutput$coefficients[1])
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))

    }

    else{

      s = fsroutput$coefficients

      mu = as.matrix(fsroutput$data$mml$mu) %*% as.numeric(s[-c(2,3)])
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))
    }

  }

  else{

    if (fsroutput$ncovs == 0) {

      mu = as.numeric(fsroutput$coefficients[2])
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))


    }

    else{

      X<-model.matrix( ~ ., data = newdata)

      s = fsroutput$coefficients

      sa = s[2]

      sb = s[-c(1,2,3)]

      sb = sb[colnames(X)]

      sb = sb[!is.na(sb)]

      sc = append(sa,sb)

      if (length(sc) != ncol(X)){

        stop('Incorrect Level Entered')


      }

      mu = as.matrix(X) %*% as.numeric(sc)
      a = exp(mu - (log(Q^(-2))*sigma)/(abs(Q)))

    }

  }



  sx = pgengamma(x,

                 mu = mu,

                 sigma = sigma,

                 Q = Q,

                 lower.tail = FALSE)

  #c = (Q^(1-2*(Q + sigma)/Q)*abs(Q))/gamma(Q^(-2))*exp(-mu)

  mx = as.numeric(a * upper_incomplete_gamma((x/a)^b, k + 1/b)/upper_incomplete_gamma((x/a)^b, k) - x)

  px = function(p){

    pc = (1-p)*sx

    px = as.numeric(qgengamma(pc,

                              mu = mu,

                              sigma = sigma,

                              Q = Q,

                              lower.tail = FALSE) - x)

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






