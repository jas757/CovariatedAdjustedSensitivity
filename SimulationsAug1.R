setwd("~/Dropbox/CovariateAdjustmentDiagnosticAccuracy/Programing/Results/")
library(rootSolve)
library(MASS)

# Unadjusted sensitivity
unad.sens = function(M,Y,A){
  sens.unad.1 = sum(Y ==1 & A ==1 & M==1)/sum(Y ==1 & A ==1)
  sens.unad.0 = mean(M[Y ==1 & A ==0])
  return(list(a.1 = sens.unad.1, a.0 = sens.unad.0))
}


# Function to calculate the model adjusted sensitivity
estimators.sens = function(data.used){
  
  # Solve the estimating equation
  mod.fit = glm(M ~.-Y, data = data.used[data.used$Y==1, ], family = "binomial")
  
  # Unadjusted estimator
  unad.est = unad.sens(data.used$M, data.used$Y, data.used$A)
  unad.est.trt = unad.est$a.1 - unad.est$a.0
  
  #Standardized estimator
  #pred.used.a.1 = exp(model.parm[1] + model.parm[2] + model.parm[3:length(model.parm)] * data.used$X)/(1 + exp(model.parm[1] + model.parm[2] + model.parm[3:length(model.parm)] * data.used$X))
  #pred.used.a.0 = exp(model.parm[1]  + model.parm[3:length(model.parm)] * data.used$X)/(1 + exp(model.parm[1]  + model.parm[3:length(model.parm)] * data.used$X))
  data.used.1 = data.used
  data.used.1$A =1 
  pred.a.1 = predict(mod.fit, newdata = data.used.1[data.used.1$Y == 1, ], type = "response")
  std.a.1 = mean(pred.a.1)
  
  # A =0
  data.used.0 = data.used
  data.used.0$A =0
  pred.a.0 = predict(mod.fit, newdata = data.used.0[data.used.0$Y == 1, ], type = "response")
  std.a.0 = mean(pred.a.0)
  
  # Influence function based estimator
  
  # Fit nuisance models
  ea = glm(A~.-Y-M, data = data.used, family = "binomial")
  ea.pred = predict(ea, newdata = data.used, type = "response")
  
  ha = glm(Y~., data = data.used, family = "binomial")
  pred.0 = data.used
  pred.0$A = 0
  pred.0$M = 1
  ha.pred.0 = predict(ha, newdata = pred.0, type = "response")
  pred.1 = data.used
  pred.1$A = 1
  pred.1$M = 1
  ha.pred.1 = predict(ha, newdata = pred.1, type = "response")
  
  ka = glm(M~.-Y, data = data.used, family = "binomial")
  pred.0 = data.used
  pred.0$A = 0
  ka.pred.0 = predict(ka, newdata = pred.0, type = "response")
  pred.1 = data.used
  pred.1$A = 1
  ka.pred.1 = predict(ka, newdata = pred.1, type = "response")
  
  la = glm(Y~.-M, data = data.used, family = "binomial")
  pred.0 = data.used
  pred.0$A = 0
  la.pred.0 = predict(la, newdata = pred.0, type = "response")
  pred.1 = data.used
  pred.1$A = 1
  la.pred.1 = predict(la, newdata = pred.1, type = "response")
  
  
  # Implement estimator
  est.inf.num.a.0 = sum(ha.pred.0 * ka.pred.0 + as.numeric(data.used$M == 1) * as.numeric(data.used$A == 0)/(1-ea.pred) * (data.used$Y - ha.pred.0) + as.numeric(data.used$A == 0)/(1-ea.pred) * ha.pred.0 * (data.used$M - ka.pred.0))
  est.inf.den.a.0 = sum(la.pred.0 + as.numeric(data.used$A == 0)/(1-ea.pred) * (data.used$Y - la.pred.0))
  est.inf.0 = est.inf.num.a.0/est.inf.den.a.0
  
  est.inf.num.a.1 = sum(ha.pred.1 * ka.pred.1 + as.numeric(data.used$M == 1) * as.numeric(data.used$A == 1)/(ea.pred) * (data.used$Y - ha.pred.1) + as.numeric(data.used$A == 1)/(ea.pred) * ha.pred.1 * (data.used$M - ka.pred.1))
  est.inf.den.a.1 = sum(la.pred.1 + as.numeric(data.used$A == 1)/(ea.pred) * (data.used$Y - la.pred.1))
  est.inf.1 = est.inf.num.a.1/est.inf.den.a.1
  
  
  return(list(model.est = mod.fit$coefficients[2], unad.est.a.1 = unad.est$a.1, unad.est.a.0 = unad.est$a.0, std.a.1 = std.a.1, std.a.0 = std.a.0, est.inf.1 = est.inf.1, est.inf.0 = est.inf.0))
}

# The estimator that requires categorical covariates
estimator.alt = function(data.used, a){
  # Fit model for P[Y=1|X, M=1, A=a]
  pred.mod.y = glm(Y~., data = data.used, family = "binomial")
  data.a.1 = data.used
  data.a.1$A = a
  data.a.1$M = 1
  pred.y = predict(pred.mod.y, type = "response", newdata = data.a.1)
  
  # Fit model for P[Y=1|X, A=a]
  pred.mod.y.no.m = glm(Y~.-M, data = data.used, family = "binomial")
  data.a.1 = data.used
  data.a.1$A = a
  pred.y.no.m = predict(pred.mod.y.no.m, type = "response", newdata = data.a.1)
  
  # Fit model for P[M=1|X, A=a]
  pred.mod.m = glm(M ~.-Y, data = data.used, family = "binomial")
  data.m.1 = data.used
  data.m.1$A = a
  pred.m = predict(pred.mod.m, type = "response", newdata = data.m.1)
  
  # Calculate alternative estimator
  num.ad = mean(pred.m * pred.y)
  den.ad = mean(pred.y.no.m)
  sens.alt = num.ad/den.ad
  
  # Calculate hybrid estimator
  sens.hyb = mean(data.used$Y[data.used$A ==a] ==1 & data.used$M[data.used$A ==a]==1)/mean(pred.y.no.m)
  return(list(sens.alt, sens.hyb))
}


# Make data for without missing data or measurement error
make.dat = function(){
  n = 500000
  A = rbinom(n,1, 0.5)
  p = 5
  covar <- matrix(0, ncol = p, nrow = p)
  rho = 0.5
  covar <- rho^abs(row(covar)-col(covar))
  mean.cov <- rep(0, p)
  xx <- mvrnorm(n, mu = mean.cov, Sigma = covar)
  X1 = xx[,1]
  X2 = xx[, 2]
  X3 = xx[, 3]
  X4 = xx[, 4]
  X5 = xx[, 5]
  prob.m =  exp(-0.1 + 0.5 * A - 1.5 *X1 + 1.5 *X2)/(1 + exp(-0.1 + 0.5 * A  - 1.5 *X1 + 1.5 *X2))
  M = rbinom(n, 1, prob.m)
  prob.y = exp(-0.1 + 0.1 * A  - 1 *X1 -1 *X2 + 1 * M)/(1 + exp(-0.1 + 0.1 * A  - 1 *X1 -1 *X2 + 1 * M))
  Y = rbinom(n, 1, prob.y)
  data.used = data.frame(X1,X2,X3,X4, X5, A,M,Y)
  data.used.bin = data.frame(X1 < mean(X1),X2 < mean(X2),X3 < mean(X3), X4 < mean(X4), X5 < mean(X5), A,M,Y)
  #summary(lm(A~X, data = data.used[data.used$Y ==1, ]))
  #Y[M==1] = rbinom(sum(M==1), 1, 0.5 + 0.1 * A + 0.1 * X)
  return(list(data.used = data.used, data.used.bin = data.used.bin))
}



# Make data for missing data setting
make.dat.miss = function(){
  n = 5000
  A = rbinom(n,1, 0.5)
  p = 5
  covar <- matrix(0, ncol = p, nrow = p)
  rho = 0.5
  covar <- rho^abs(row(covar)-col(covar))
  mean.cov <- rep(0, p)
  xx <- mvrnorm(n, mu = mean.cov, Sigma = covar)
  X1 = xx[,1]
  X2 = xx[, 2]
  X3 = xx[, 3]
  X4 = xx[, 4]
  X5 = xx[, 5]
  prob.m =  exp(-0.1 + 0.5 * A - 1.5 *X1 + 1.5 *X2)/(1 + exp(-0.1 + 0.5 * A  - 1.5 *X1 + 1.5 *X2))
  M = rbinom(n, 1, prob.m)
  prob.y = exp(-0.1 + 0.1 * A  - 1 *X1 -1 *X2 + 1 * M)/(1 + exp(-0.1 + 0.1 * A  - 1 *X1 -1 *X2 + 1 * M))
  Y = rbinom(n, 1, prob.y)
  R = rbinom(n, 1, 0.1 + 0.3 * as.numeric(X1 > 0.5))
  X1[R==1] = mean(X1[R==0])
  R = rbinom(n, 1, 0.1 + 0.3 * as.numeric(X2 > 0.5))
  X2[R==1] = mean(X2[R==0])
  R = rbinom(n, 1, 0.1 + 0.3 * as.numeric(X3 > 0.5))
  X3[R==1] = mean(X3[R==0])
  R = rbinom(n, 1, 0.1 + 0.3 * as.numeric(X4 > 0.5))
  X4[R==1] = mean(X4[R==0])
  R = rbinom(n, 1, 0.1 + 0.3 * as.numeric(X5 > 0.5))
  X5[R==1] = mean(X5[R==0])

  data.used = data.frame(X1,X2,X3,X4, X5, A,M,Y)
  data.used.bin = data.frame(X1 < mean(X1),X2 < mean(X2),X3 < mean(X3), X4 < mean(X4), X5 < mean(X5), A,M,Y)
  #summary(lm(A~X, data = data.used[data.used$Y ==1, ]))
  #Y[M==1] = rbinom(sum(M==1), 1, 0.5 + 0.1 * A + 0.1 * X)
  return(list(data.used = data.used, data.used.bin = data.used.bin))
}




# Make data for measurement error
make.dat.meas.err = function(){
    n = 5000
    A = rbinom(n,1, 0.5)
    p = 5
    covar <- matrix(0, ncol = p, nrow = p)
    rho = 0.5
    covar <- rho^abs(row(covar)-col(covar))
    mean.cov <- rep(0, p)
    xx <- mvrnorm(n, mu = mean.cov, Sigma = covar)
    X1 = xx[,1]
    X2 = xx[, 2]
    X3 = xx[, 3]
    X4 = xx[, 4]
    X5 = xx[, 5]
    prob.m =  exp(-0.1 + 0.5 * A - 1.5 *X1 + 1.5 *X2)/(1 + exp(-0.1 + 0.5 * A  - 1.5 *X1 + 1.5 *X2))
    M = rbinom(n, 1, prob.m)
    prob.y = exp(-0.1 + 0.1 * A  - 1 *X1 -1 *X2 + 1 * M)/(1 + exp(-0.1 + 0.1 * A  - 1 *X1 -1 *X2 + 1 * M))
    Y = rbinom(n, 1, prob.y)
    X1 = X1 + rnorm(n, 0.1, 0.1)
    X2 = X2 + rnorm(n, 0.1, 0.1)
    X3 = X3 + rnorm(n, 0.1, 0.1)
    X4 = X4 + rnorm(n, 0.1, 0.1)
    X5 = X5 + rnorm(n, 0.1, 0.1)
    data.used = data.frame(X1,X2,X3,X4, X5, A,M,Y)
    data.used.bin = data.frame(X1 < mean(X1),X2 < mean(X2),X3 < mean(X3), X4 < mean(X4), X5 < mean(X5), A,M,Y)
    #summary(lm(A~X, data = data.used[data.used$Y ==1, ]))
    #Y[M==1] = rbinom(sum(M==1), 1, 0.5 + 0.1 * A + 0.1 * X)
    return(list(data.used = data.used, data.used.bin = data.used.bin))
}



# Simulate data
n.sim = 1000
sens.unad.1 = rep(NA, n.sim)
sens.unad.0 = rep(NA, n.sim)
sens.alt.1 = rep(NA, n.sim)
sens.alt.0 = rep(NA, n.sim)
sens.hyb.1 = rep(NA, n.sim)
sens.hyb.0 = rep(NA, n.sim)
std.1 = rep(NA, n.sim)
std.0 = rep(NA, n.sim)
inf.1 = rep(NA, n.sim)
inf.0 = rep(NA, n.sim)
mod.est = rep(NA, n.sim)
set.seed = 1
for(i in 1:n.sim){
  data.used = make.dat()
  
  # Fit the estimators
  est.sens = estimators.sens(data.used$data.used)
  sens.unad.1[i] = est.sens$unad.est.a.1
  sens.unad.0[i] = est.sens$unad.est.a.0
  std.1[i] = est.sens$std.a.1
  std.0[i] = est.sens$std.a.0
  mod.est[i] = est.sens$model.est
  
  # Alternative estimator
  est.sens.alt.1 = estimator.alt(data.used$data.used.bin, a=1)
  est.sens.alt.0 = estimator.alt(data.used$data.used.bin, a=0)
  sens.alt.1[i] = est.sens.alt.1[[1]]
  sens.alt.0[i] = est.sens.alt.0[[1]]
  sens.hyb.1[i] = est.sens.alt.1[[2]]
  sens.hyb.0[i] = est.sens.alt.0[[2]]
  
  # Influence based
  inf.0[i] = est.sens$est.inf.0
  inf.1[i] = est.sens$est.inf.1
}

data.comb = list(sens.unad.1 = sens.unad.1,
                 sens.unad.0 = sens.unad.0,
                 sens.alt.1 = sens.alt.1,
                 sens.alt.0 = sens.alt.0,
                 sens.hyb.1 = sens.hyb.1,
                 sens.hyb.0 = sens.hyb.0,
                 std.1 = std.1,
                 std.0 = std.0,
                 inf.1 = inf.1,
                 inf.0 = inf.0)

save(data.comb, file = "ResultsNoMissMeasn5000.rda")

# Influence vs unadjusted
mean(inf.1) - mean(sens.unad.1)
mean(inf.0) - mean(sens.unad.0)
var(sens.unad.0)/var(inf.0)
var(sens.unad.1)/var(inf.1)


# Adjusted vs unadjusted
mean(std.1) - mean(sens.unad.1)
mean(std.0) - mean(sens.unad.0)
var(sens.unad.0)/var(sens.alt.0)
var(sens.unad.1)/var(sens.alt.1)

# Hybrid vs unadjusted
mean(sens.hyb.1) - mean(sens.unad.1)
mean(sens.hyb.0) - mean(sens.unad.0)
var(sens.unad.0)/var(sens.hyb.0)
var(sens.unad.0)/var(std.0)

# Standardized vs unadjusted
mean(std.1) - mean(sens.unad.1)
mean(std.0) - mean(sens.unad.0)
var(sens.unad.1)/var(std.1)
var(sens.unad.0)/var(std.0)
