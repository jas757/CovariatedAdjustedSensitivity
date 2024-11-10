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
  

  return(list(sens.alt))
}

