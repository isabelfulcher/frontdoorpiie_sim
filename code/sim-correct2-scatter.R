############################################################
# SIMULATIONS -- CLUSTER -- CORRECT (LEAVE OUT CONFOUNDER) #
############################################################

library(readr)
library(foreach)
library(MASS)
library(doParallel)
library(magrittr)
library(numDeriv)

cl <- makeCluster(10) 
registerDoParallel(cl)

source("code/Simulation-continuous-functions_v1.R")
load("truth/frontdoor-truth.rda")


set.seed(400)

S=10000
n=1000

parOut <- foreach(s=1:S) %dopar% {
  data <- gen.med.data.continuous(n,p.c1,p.c3,alpha.truth,omega.truth,beta.truth,theta.truth,sigma.m,sigma.y)
  
  fit.z <- lm(m ~ a + c1 + c2 + I(c1*c2),data=data) 
  fit.y <- lm(y ~ a + m + I(a*m) + c1 + c2 + I(c1*c2),data=data) 
  fit.a <- glm(a ~ c1 + c2 + I(c1*c2),data=data,family=binomial)
  
  alpha.hat <- summary(fit.a)$coefficients[,1]
  beta.hat <- summary(fit.z)$coefficients[,1]
  theta.hat <- summary(fit.y)$coefficients[,1]
  
  confounders <- cbind(data[,4:5],data[,4]*data[,5])
  
  # MLE #
  out.mle <- piie.mle.variance.function.cont(confounders,data$a,theta.hat,beta.hat,alpha.hat,0,1)
  prop.bias.mle <- (out.mle[2]-truth.est[3])/truth.est[3]
  
  ci.low <- out.mle[2] - qnorm(.975)*sqrt(out.mle[3])
  ci.up <- out.mle[2] + qnorm(.975)*sqrt(out.mle[3])
  ci.coverage.mle <- as.numeric(truth.est[3] >= ci.low & truth.est[3] <= ci.up)
  
  out.mle <- cbind(out.mle,prop.bias.mle,ci.coverage.mle)
  
  # SP 1 #
  out.sp.1 <- piie.sp.1.variance.function.cont(data$m,data$y,fit.z,0,1)
  prop.bias.sp.1 <- (out.sp.1[2]-truth.est[3])/truth.est[3]
  
  ci.low <- out.sp.1[2] - qnorm(.975)*sqrt(out.sp.1[3])
  ci.up <- out.sp.1[2] + qnorm(.975)*sqrt(out.sp.1[3])
  ci.coverage.sp.1 <- as.numeric(truth.est[3] >= ci.low & truth.est[3] <= ci.up)
  
  out.sp.1 <- cbind(out.sp.1,prop.bias.sp.1,ci.coverage.sp.1)
  
  # SP 2 #
  out.sp.2 <- piie.sp.2.variance.function.cont(confounders,data$a,data$y,c(1,1,1),fit.a,fit.y,0,1)
  prop.bias.sp.2 <- (out.sp.2[2]-truth.est[3])/truth.est[3]
  
  ci.low <- out.sp.2[2] - qnorm(.975)*sqrt(out.sp.2[3])
  ci.up <- out.sp.2[2] + qnorm(.975)*sqrt(out.sp.2[3])
  ci.coverage.sp.2 <- as.numeric(truth.est[3] >= ci.low & truth.est[3] <= ci.up)
  
  out.sp.2 <- cbind(out.sp.2,prop.bias.sp.2,ci.coverage.sp.2)
  
  # DR #
  out.sp <- piie.sp.variance.function.cont(confounders,data$a,data$m,data$y,c(1,1,1),c(1,1,1),c(1,1,1),fit.a,fit.z,fit.y,0,1)
  prop.bias.sp <- (out.sp[2]-truth.est[3])/truth.est[3]
  
  ci.low <- out.sp[2] - qnorm(.975)*sqrt(out.sp[3])
  ci.up <- out.sp[2] + qnorm(.975)*sqrt(out.sp[3])
  ci.coverage.sp <- as.numeric(truth.est[3] >= ci.low & truth.est[3] <= ci.up)
  
  out.sp <- cbind(out.sp,prop.bias.sp,ci.coverage.sp)
  
  mean.y <- mean(data$y)
  
  list(out.mle,out.sp.1,out.sp.2,out.sp,mean.y)
  
}

getter <- function(iter, index) iter[[index]]

mle <- sapply(parOut, getter, 1) %>% t()
sp.1 <- sapply(parOut, getter, 2) %>% t()
sp.2 <- sapply(parOut, getter, 3) %>% t()
sp <- sapply(parOut, getter, 4) %>% t()
mean.y <- sapply(parOut, getter, 5) %>% t()


outlist <- list(mle,sp.1,sp.2,sp,mean.y)
saveRDS(outlist, file = paste0("output/out_correct2.rds"))
