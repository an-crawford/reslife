###Written by Ka Lok Lee

###Load flexsurv
library(flexsurv)

#data(package='survival')

###Check dataset bc
head(bc)

###Run Weibull
fsw<-flexsurvreg(formula = Surv(recyrs, censrec) ~ group, data = bc,dist = "weibull")
fsw$coefficients #the other way
fsw$res #positive-enforced parameters

fsw$ncovs #number of covariates

fsw$dlist$name

###Run Exponential
fse<-flexsurvreg(formula = Surv(recyrs, censrec) ~ 1, data = bc,dist = "exp")
fse$dlist$name
fse$coefficients


###Exponential MRL function
exp_mrl <- function(fsoutput,life) {
  
  if (fsoutput$dlist$name != 'exp') {
    print("error")
  }
  
  else  {
    
    if (fsoutput$ncovs == 0) {
      
      mean_res_life = 1/exp(fsoutput$coefficients)
      
    }
    
    else {
      mean_res_life = 1/exp(as.matrix(fsoutput$data$mml$rate) %*% as.numeric(fsoutput$coefficients));
    }
    return (mean_res_life)
    
  }
}

###Calculate MRL
unique(exp_mrl(fse,1))

### plot with multiple "life's"

life_vector = seq(1:10)
exp_mrl_many = seq(1:10)
for (i in life_vector) {
  
  exp_mrl_many[i] = exp_mrl(fse,i)
  
}


plot(life_vector, exp_mrl_many, type="l")

###Calculate using the wrong object
exp_mrl(fsw,1)



###Simulation

N=100000

b0 <- 0.5
b1 <- 1
b2 <- -.2

x0 <- rep(1,N)
x1 <- rnorm(N, .1, .2)
x2 <- rbinom(N,1,.6)

rr <- exp(cbind(x0,x1,x2) %*% c(b0,b1,b2))
length(rr)
time <- rexp(N, rate=rr)
length(time)
mean(time)

time_2 = pmin(.5,time)
event = ifelse(time_2<.5,1,0)

exp_df <- data.frame( t = time_2,
                  status = event,
                  x1 = x1,
                  x2 = x2)

head(exp_df)

fs_s<-flexsurvreg(formula = Surv(t, status) ~ x1 +x2 , data = exp_df,dist = "exp")

fs_s$coefficients




fs_s<-flexsurvreg(formula = Surv(t, status) ~ 1 , data = exp_df,dist = "exp")
fs_s$ncovs
fs_s$res
fs_s$coefficients


