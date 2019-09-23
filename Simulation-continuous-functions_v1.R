#################################
### FUNCITIONS FOR SIMULATION ###
#################################

expit <- function(x){
  output <- exp(x) / (1 + exp(x))
  return(output)
}

######## GENERATE DATASET ###########
gen.med.data.continuous <- function(n,p.c1,p.c3,alpha,omega,beta,theta,sigma.m,sigma.y){
  
  #Generate first confounder 
  c1 <- rbinom(n,1,p.c1)
  
  #Generate second confounder (given c1)
  p.c2 <- expit(cbind(rep(1,n),c1)%*%alpha)*(1-as.numeric(alpha[1] == 0 &  alpha[2] == 0)) 
  c2 <- rbinom(n,1,p.c2) 
  
  #Generate a confounder of A & Y (independent of other confounders)
  c3 <- rbinom(n,1,p.c3)
  
  #Generate binary exposure (given c1, c2 and c1*c2)
  p.a <- expit(cbind(rep(1,n),c1,c2,c1*c2,c3)%*%omega)
  a <- rbinom(n,1,p.a) 
  
  #Generate continuous mediator (given a, c1, c2 and all possible combinations)
  mean.m <- cbind(rep(1,n),a,c1,c2,c1*c2,c3)%*%beta
  m <- rnorm(n,mean.m,sigma.m)
  
  #Generate continuous outcome Y (given a,m,c1,c2, c3 and possible combos)
  mean.y <- cbind(rep(1,n),a,m,a*m,c1,c2,c1*c2,c3)%*%theta
  y <- rnorm(n,mean.y,sigma.y)
  
  sim.data <- data.frame(cbind(y,a,m,c1,c2,c3))
  
  return(sim.data)
  
}

####### CALCULATE TRUE VALUE OF PIIE AND MEAN ####### 

psi.truth.function.continuous <- function(n,p.c1,p.c3,alpha,omega,beta.truth,theta.truth){
  
  p.c2 <- expit(sum(alpha))*p.c1 + expit(alpha[1])*(1-p.c1) #E(C2)
  
  p.a <-  ( ((expit(omega[1]+omega[5])*(1-p.c1)*(1-expit(alpha[1]))
              + expit(sum(omega[1:2]+omega[5]))*p.c1*(1-expit(sum(alpha)))
              + expit(omega[1] + omega[3]+omega[5])*(1-p.c1)*expit(alpha[1])
              + expit(sum(omega))*p.c1*expit(sum(alpha))))*p.c3 
            + ((expit(omega[1])*(1-p.c1)*(1-expit(alpha[1]))
                + expit(sum(omega[1:2]))*p.c1*(1-expit(sum(alpha)))
                + expit(omega[1] + omega[3])*(1-p.c1)*expit(alpha[1])
                + expit(sum(omega[1:4]))*p.c1*expit(sum(alpha))))*(1-p.c3) ) #E(A)
  
  p.c1c2 <- p.c1*expit(sum(alpha)) #E(C1*C2) 
  
  p.ac1 <- p.c1*( expit(sum(omega))*p.c3*expit(sum(alpha))
                  + expit(sum(omega[1:4]))*(1-p.c3)*expit(sum(alpha))
                  + expit(omega[1]+omega[2]+omega[5])*p.c3*(1-expit(sum(alpha)))
                  + expit(omega[1]+omega[2])*(1-p.c3)*(1-expit(sum(alpha))) ) #E(AC1)
  p.ac2 <- ( expit(sum(omega))*p.c3*p.c1c2
             + expit(sum(omega[1:4]))*(1-p.c3)*p.c1c2
             + expit(omega[1] + omega[3] + omega[5])*p.c3*(1-p.c1)*expit(alpha[1])
             + expit(omega[1] + omega[3])*(1-p.c3)*(1-p.c1)*expit(alpha[1]) ) #E(AC2) - something is up
  
  p.ac3 <- p.c3*( expit(sum(omega))*p.c1*p.c2
                  + expit(omega[1] + omega[2] + omega[5])*p.c1*(1-p.c2)
                  + expit(omega[1]+omega[3]+omega[5])*(1-p.c1)*p.c2
                  + expit(omega[1]+omega[2])*(1-p.c1)*(1-p.c2)) #E(A*C3)
  
  p.ac1c2 <- p.c1c2*(expit(sum(omega))*p.c3 + 
                       expit(sum(omega[1:4]))*(1-p.c3)) #E(A*C1*C2)
  
  psi <- (theta.truth[1] + theta.truth[3]*beta.truth[1] 
                + (theta.truth[2] + theta.truth[4]*beta.truth[1])*p.a
                + (theta.truth[3]*beta.truth[3:length(beta.truth)] + theta.truth[5:length(theta.truth)])%*%rbind(p.c1,p.c2,p.c1c2,p.c3)
                + theta.truth[4]*beta.truth[3:length(beta.truth)]%*%rbind(p.ac1,p.ac2,p.ac1c2,p.ac3))
  
  mean.y <- (theta.truth[1] + theta.truth[3]*beta.truth[1] 
             + (theta.truth[2] + theta.truth[3]*beta.truth[2] + theta.truth[4]*beta.truth[1])*p.a 
             + (theta.truth[3]*beta.truth[3:length(beta.truth)] + theta.truth[5:length(theta.truth)])%*%rbind(p.c1,p.c2,p.c1c2,p.c3)
             + theta.truth[4]*beta.truth[3:length(beta.truth)]%*%rbind(p.ac1,p.ac2,p.ac1c2,p.ac3)
             + theta.truth[4]*beta.truth[2]*p.a)
  
  piie.1 <- mean.y - psi
  
  piie.2 <- theta.truth[3]*beta.truth[2]*p.a + theta.truth[4]*beta.truth[2]*p.a
  
  output <- matrix(c(mean.y,psi,piie.1,piie.2),1,4)
  colnames(output) <- c("Mean Y","PSI","PIIE - diff","PIIE - exp")
  
  return(output)
}


###### FUNCTION THAT CALCULATE PSI MLE###### 
psi.mle.function.cont <- function(cov.vals.all,exposure.data,beta.hat,theta.hat,astar,interaction){
  
  mean.cov.all <- apply(cov.vals.all,2,mean)
  mean.exposure <- mean(exposure.data)
  mean.cov.exposure <- apply(exposure.data*cov.vals.all,2,mean)
  
  
  if (interaction == 1){psi <- (theta.hat[1] + theta.hat[3]*beta.hat[1] + theta.hat[3]*beta.hat[2]*astar
                                + (theta.hat[2] + theta.hat[4]*beta.hat[1] + theta.hat[4]*beta.hat[2]*astar)*mean.exposure
                                + (theta.hat[3]*beta.hat[3:length(beta.hat)] + theta.hat[5:length(theta.hat)])%*%t(t(mean.cov.all))
                                + theta.hat[4]*beta.hat[3:length(beta.hat)]%*%t(t(mean.cov.exposure)))
                } else {psi <- (theta.hat[1] + theta.hat[3]*beta.hat[1] + theta.hat[3]*beta.hat[2]*astar
                                + theta.hat[2]*mean.exposure
                                + (theta.hat[3]*beta.hat[3:length(beta.hat)] + theta.hat[4:length(theta.hat)])%*%t(t(mean.cov.all)))}

  return(psi)
}

piie.mle.variance.function.cont <- function(cov.vals.all,exposure.data,theta.hat,beta.hat,alpha.hat,astar,interaction){
  
  n <- length(exposure.data)
  
  mean.cov.all <- apply(cov.vals.all,2,mean)
  mean.exposure <- mean(exposure.data)
  mean.cov.exposure <- apply(exposure.data*cov.vals.all,2,mean)
  
  
  if (interaction == 1){psi <- (theta.hat[1] + theta.hat[3]*beta.hat[1] + theta.hat[3]*beta.hat[2]*astar
                                + (theta.hat[2] + theta.hat[4]*beta.hat[1] + theta.hat[4]*beta.hat[2]*astar)*mean.exposure
                                + (theta.hat[3]*beta.hat[3:length(beta.hat)] + theta.hat[5:length(theta.hat)])%*%t(t(mean.cov.all))
                                + theta.hat[4]*beta.hat[3:length(beta.hat)]%*%t(t(mean.cov.exposure)))
          } else {psi <- (theta.hat[1] + theta.hat[3]*beta.hat[1] + theta.hat[3]*beta.hat[2]*astar
                          + theta.hat[2]*mean.exposure
                          + (theta.hat[3]*beta.hat[3:length(beta.hat)] + theta.hat[4:length(theta.hat)])%*%t(t(mean.cov.all)))}
  
  if(interaction == 1){est.piie <- theta.hat[3]*beta.hat[2]*mean.exposure + theta.hat[4]*beta.hat[2]*mean.exposure
              } else {est.piie <- theta.hat[3]*beta.hat[2]*mean.exposure   }
  
  if(interaction ==1){var.piie <- ( (var(exposure.data)/n)*(beta.hat[2]^2)*(theta.hat[3] + theta.hat[4])^2 
                  + (mean.exposure^2)*((theta.hat[3] + theta.hat[4])^2)*vcov(fit.z)[2,2]
                  + mean.exposure*beta.hat[2]*(mean.exposure*beta.hat[2]*vcov(fit.y)[3,3] + mean.exposure*beta.hat[2]*vcov(fit.y)[3,4])
                  + mean.exposure*beta.hat[2]*(mean.exposure*beta.hat[2]*vcov(fit.y)[3,4] + mean.exposure*beta.hat[2]*vcov(fit.y)[4,4]) )
              } else { var.piie <- ( ((mean.exposure*theta.hat[3])^2)*vcov(fit.z)[2,2]
                                    + ((mean.exposure*beta.hat[2])^2)*vcov(fit.y)[3,3]
                                    + ((beta.hat[2]*theta.hat[3])^2)*var(exposure.data)/n )  }
  
  output <- cbind(psi,est.piie,var.piie)
  colnames(output) <- c("PSI","PIIE","Var PIIE")
  
  return(output)
  
}


###### FUNCTION THAT CALCULATE PSI DOUBLY ROBUST SP###### 
psi.sp.function.cont <- function(cov.vals.all,outcome,i.y,i.z,i.a,fit.a,fit.z,fit.y,astar,interaction){
 
  n <- length(outcome)
  
  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))
  
  sigma <- summary(fit.z)$sigma
  
  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[,2] <- 0
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  theta.hat <- summary(fit.y)$coefficients[,1]
  beta.hat <- summary(fit.z)$coefficients[,1]
  alpha.hat <- summary(fit.a)$coefficients[,1]
  
  cov.vals.y <- cov.vals.all[,which(i.y==1)]
  cov.vals.z <- cov.vals.all[,which(i.z==1)]
  cov.vals.a <- cov.vals.all[,which(i.a==1)]
 
  z.mean_astar <- model.matrix.z_astar%*%beta.hat
  z.mean_ind <- model.matrix.z%*%beta.hat
  a.mean <- expit(model.matrix.a%*%alpha.hat)
  y.mean <- model.matrix.y%*%theta.hat

  if (interaction == 1){sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],a.mean*model.matrix.y[,3],cov.vals.y))%*%theta.hat
  } else { sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],cov.vals.y))%*%theta.hat }
  
  if (interaction == 1){sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,2],z.mean_astar,model.matrix.y[,2]*z.mean_astar,cov.vals.y))%*%theta.hat
  } else { sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,2],z.mean_astar,cov.vals.y))%*%theta.hat}
  
  if (interaction == 1){sum.az <- as.matrix(cbind(rep(1,n),a.mean,z.mean_ind,a.mean*z.mean_ind,cov.vals.y))%*%theta.hat
  } else {sum.az <- as.matrix(cbind(rep(1,n),a.mean,z.mean_ind,cov.vals.y))%*%theta.hat}


  psi <- sum( (outcome - y.mean)*
                   (dnorm(model.matrix.y[,3],z.mean_astar,sigma)/dnorm(model.matrix.y[,3],z.mean_ind,sigma))
                 + ((1-model.matrix.y[,2])/(1-a.mean))*(sum.a - sum.az) 
                 + sum.z )/nrow(model.matrix.y)
  
  return(psi)

}

piie.sp.variance.function.cont <- function(cov.vals.all,exposure,intermediate,outcome,i.y,i.z,i.a,fit.a,fit.z,fit.y,astar,interaction){
 
  n <- length(exposure)
  
  sigma <- summary(fit.z)$sigma
  
  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))
  
  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[,2] <- 0
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  theta.hat <- summary(fit.y)$coefficients[,1]
  beta.hat <- summary(fit.z)$coefficients[,1]
  alpha.hat <- summary(fit.a)$coefficients[,1]
  
  cov.vals.y <- cov.vals.all[,which(i.y==1)]
  cov.vals.z <- cov.vals.all[,which(i.z==1)]
  cov.vals.a <- cov.vals.all[,which(i.a==1)]
  
  z.mean_astar <- model.matrix.z_astar%*%beta.hat
  z.mean_ind <- model.matrix.z%*%beta.hat
  a.mean <- expit(model.matrix.a%*%alpha.hat)
  y.mean <- model.matrix.y%*%theta.hat
  
  if (interaction == 1){sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],a.mean*model.matrix.y[,3],cov.vals.y))%*%theta.hat
  } else { sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],cov.vals.y))%*%theta.hat }
  
  if (interaction == 1){sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,2],z.mean_astar,model.matrix.y[,2]*z.mean_astar,cov.vals.y))%*%theta.hat
  } else { sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,2],z.mean_astar,cov.vals.y))%*%theta.hat}
  
  if (interaction == 1){sum.az <- as.matrix(cbind(rep(1,n),a.mean,z.mean_ind,a.mean*z.mean_ind,cov.vals.y))%*%theta.hat
  } else {sum.az <- as.matrix(cbind(rep(1,n),a.mean,z.mean_ind,cov.vals.y))%*%theta.hat}
  
  psi.sp.ind <-  ((outcome - y.mean)*
                  (dnorm(model.matrix.y[,3],z.mean_astar,sigma)/dnorm(model.matrix.y[,3],z.mean_ind,sigma))
                  + ((1-model.matrix.y[,2])/(1-a.mean))*(sum.a - sum.az) 
                  + sum.z )
  
  piie.sp.ind <- outcome - psi.sp.ind
  
  piie.sp <- mean(outcome) - mean(psi.sp.ind)
  
  score.sp <- cbind( model.matrix.a*c(exposure - a.mean),
                     model.matrix.z*c(intermediate - z.mean_ind),
                     model.matrix.y*c(outcome - y.mean),
                     (piie.sp.ind - piie.sp))
  
  estimates <- c(alpha.hat,beta.hat,theta.hat)
  len.a <- length(coefficients(fit.a))
  len.z <- length(coefficients(fit.z))
  
  #take the derivative and plugs in inputs
  deriv.sp <- numDeriv::jacobian(U.sp,c(estimates,piie.sp),model.matrix.a=model.matrix.a,model.matrix.z=model.matrix.z,model.matrix.y=model.matrix.y,data.y=outcome,data.z=intermediate,data.a=exposure,i.y=i.y,i.z=i.z,i.a=i.a,cov.vals.all=cov.vals.all,len.a=len.a,len.z=len.z,n=n,sigma=sigma,interaction=interaction)
  
  #Calculate variance matrix
  var.sp <- (solve(deriv.sp)%*%t(score.sp)%*%score.sp%*%t(solve(deriv.sp)))
  
  #Variance for our estimator 
  piie.var.sp <- var.sp[length(c(estimates,piie.sp)),length(c(estimates,piie.sp))]
  
  output <- cbind(mean(psi.sp.ind),piie.sp,piie.var.sp)
  colnames(output) <- c("PSI","PIIE","Var PIIE")
  
  return(output)
  
}
U.sp <- function(estimates,model.matrix.a,model.matrix.z,model.matrix.y,data.a,data.z,data.y,i.y,i.z,i.a,cov.vals.all,len.a,len.z,n,sigma,interaction){
  
  alpha.hat <- estimates[1:len.a]
  beta.hat <- estimates[(len.a+1):(len.a+len.z)]
  theta.hat <- estimates[(len.a+len.z+1):(length(estimates) - 1)]
  piie.est <- estimates[length(estimates)]
  
  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[,2] <- 0
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  cov.vals.y <- cov.vals.all[,which(i.y==1)]
  cov.vals.z <- cov.vals.all[,which(i.z==1)]
  cov.vals.a <- cov.vals.all[,which(i.a==1)]
  
  z.mean_astar <- model.matrix.z_astar%*%beta.hat
  z.mean_ind <- model.matrix.z%*%beta.hat
  a.mean <- expit(model.matrix.a%*%alpha.hat)
  y.mean <- model.matrix.y%*%theta.hat
  
  if (interaction == 1){sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],a.mean*model.matrix.y[,3],cov.vals.y))%*%theta.hat
  } else { sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],cov.vals.y))%*%theta.hat }
  
  if (interaction == 1){sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,2],z.mean_astar,model.matrix.y[,2]*z.mean_astar,cov.vals.y))%*%theta.hat
  } else { sum.z <- as.matrix(cbind(rep(1,n),model.matrix.y[,2],z.mean_astar,cov.vals.y))%*%theta.hat }
  
  if (interaction == 1){sum.az <- as.matrix(cbind(rep(1,n),a.mean,z.mean_ind,a.mean*z.mean_ind,cov.vals.y))%*%theta.hat
  } else {sum.az <- as.matrix(cbind(rep(1,n),a.mean,z.mean_ind,cov.vals.y))%*%theta.hat}
  
  psi.sp.ind <- ( (data.y - y.mean)*
                        (dnorm(model.matrix.y[,3],z.mean_astar,sigma)/dnorm(model.matrix.y[,3],z.mean_ind,sigma))
                        + ((1-model.matrix.y[,2])/(1-a.mean))*(sum.a - sum.az) 
                        + sum.z )
  
  piie.sp.ind <- data.y - psi.sp.ind
  
  piie.sp <- mean(data.y) - mean(psi.sp.ind)
  
  score.sp <- cbind( model.matrix.a*c(data.a - a.mean),
                     model.matrix.z*c(data.z - z.mean_ind),
                     model.matrix.y*c(data.y - y.mean),
                     (piie.sp.ind - piie.est))
  
  deriv <- matrix(1,1,n)%*%score.sp
  
  return(deriv)
}
  
###### FUNCTION THAT CALCULATES PSI SP 1 (Model for Z) ######

piie.sp.1.variance.function.cont <- function(intermediate,outcome,fit.z,astar,interaction){
  
  n <- length(outcome)
  
  sigma <- summary(fit.z)$sigma
  
  model.matrix.z <- as.matrix(data.frame(model.matrix(fit.z)))

  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[,2] <- 0
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  beta.hat <- summary(fit.z)$coefficients[,1]

  z.mean_astar <- model.matrix.z_astar%*%beta.hat
  z.mean_ind <- model.matrix.z%*%beta.hat

  psi.sp.ind <-  outcome*(dnorm(intermediate,z.mean_astar,sigma)/dnorm(intermediate,z.mean_ind,sigma))
  
  piie.sp.ind <- outcome - psi.sp.ind
  
  piie.sp <- mean(outcome) - mean(psi.sp.ind)
  
  score.sp <- cbind(model.matrix.z*c(intermediate - z.mean_ind),
                     (piie.sp.ind - piie.sp))
  
  estimates <- beta.hat
  len.z <- length(coefficients(fit.z))
  
  #take the derivative and plugs in inputs
  deriv.sp <- numDeriv::jacobian(U.sp.1,c(estimates,piie.sp),model.matrix.z=model.matrix.z,data.y=outcome,data.z=intermediate,len.z=len.z,n=n,sigma=sigma,interaction=interaction)
  
  #Calculate variance matrix
  var.sp <- (solve(deriv.sp)%*%t(score.sp)%*%score.sp%*%t(solve(deriv.sp)))
  
  #Variance for our estimator 
  piie.var.sp <- var.sp[length(c(estimates,piie.sp)),length(c(estimates,piie.sp))]
  
  output <- cbind(mean(psi.sp.ind),piie.sp,piie.var.sp)
  colnames(output) <- c("PSI.1","PIIE.1","Var PIIE.1")
  
  return(output)
  
}
U.sp.1 <- function(estimates,model.matrix.z,data.y,data.z,len.z,n,sigma,interaction){
  
  beta.hat <- estimates[1:len.z]
  piie.est <- estimates[length(estimates)]
  
  model.matrix.z_astar <- data.frame(model.matrix.z)
  model.matrix.z_astar[,2] <- 0
  model.matrix.z_astar <- as.matrix(model.matrix.z_astar)
  
  z.mean_astar <- model.matrix.z_astar%*%beta.hat
  z.mean_ind <- model.matrix.z%*%beta.hat
  
  psi.sp.ind <- data.y*(dnorm(data.z,z.mean_astar,sigma)/dnorm(data.z,z.mean_ind,sigma))
  
  piie.sp.ind <- data.y - psi.sp.ind
  
  piie.sp <- mean(data.y) - mean(psi.sp.ind)
  
  score.sp <- cbind(model.matrix.z*c(data.z - z.mean_ind),
                     (piie.sp.ind - piie.est))
  
  deriv <- matrix(1,1,n)%*%score.sp
  
  return(deriv)
}



###### FUNCTION THAT CALCULATES PSI SP 2 (MODEL FOR A,Y) ######

piie.sp.2.variance.function.cont <- function(cov.vals.all,exposure,outcome,i.y,fit.a,fit.y,astar,interaction){
  
  n <- length(exposure)
  
  model.matrix.a <- as.matrix(data.frame(model.matrix(fit.a)))
  model.matrix.y <- as.matrix(data.frame(model.matrix(fit.y)))

  theta.hat <- summary(fit.y)$coefficients[,1]
  alpha.hat <- summary(fit.a)$coefficients[,1]
  
  cov.vals.y <- cov.vals.all[,which(i.y==1)]
  
  a.mean <- expit(model.matrix.a%*%alpha.hat)
  y.mean <- model.matrix.y%*%theta.hat
  
  
  if (interaction == 1){sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],a.mean*model.matrix.y[,3],cov.vals.y))%*%theta.hat
  } else { sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],cov.vals.y))%*%theta.hat }
  
  psi.sp.ind <-  ((1-model.matrix.y[,2])/(1-a.mean))*sum.a 
  
  piie.sp.ind <- outcome - psi.sp.ind
  
  piie.sp <- mean(outcome) - mean(psi.sp.ind)
  
  score.sp <- cbind( model.matrix.a*c(exposure - a.mean),
                     model.matrix.y*c(outcome - y.mean),
                     (piie.sp.ind - piie.sp))
  
  estimates <- c(alpha.hat,theta.hat)
  len.a <- length(coefficients(fit.a))

  #take the derivative and plugs in inputs
  deriv.sp <- numDeriv::jacobian(U.sp.2,c(estimates,piie.sp),model.matrix.a=model.matrix.a,model.matrix.y=model.matrix.y,data.y=outcome,data.a=exposure,i.y=i.y,cov.vals.all=cov.vals.all,len.a=len.a,n=n,sigma=sigma,interaction=interaction)
  
  #Calculate variance matrix
  var.sp <- (solve(deriv.sp)%*%t(score.sp)%*%score.sp%*%t(solve(deriv.sp)))
  
  #Variance for our estimator 
  piie.var.sp <- var.sp[length(c(estimates,piie.sp)),length(c(estimates,piie.sp))]
  
  output <- cbind(mean(psi.sp.ind),piie.sp,piie.var.sp)
  colnames(output) <- c("PSI.2","PIIE.2","Var PIIE.2")
  
  return(output)
  
}
U.sp.2 <- function(estimates,model.matrix.a,model.matrix.y,data.a,data.y,i.y,cov.vals.all,len.a,n,sigma,interaction){
  
  alpha.hat <- estimates[1:len.a]
  theta.hat <- estimates[(len.a+1):(length(estimates) - 1)]
  piie.est <- estimates[length(estimates)]
  
  cov.vals.y <- cov.vals.all[,which(i.y==1)]

  a.mean <- expit(model.matrix.a%*%alpha.hat)
  y.mean <- model.matrix.y%*%theta.hat
  
  if (interaction == 1){sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],a.mean*model.matrix.y[,3],cov.vals.y))%*%theta.hat
  } else { sum.a <- as.matrix(cbind(rep(1,n),a.mean,model.matrix.y[,3],cov.vals.y))%*%theta.hat }
  
  psi.sp.ind <- ((1-model.matrix.y[,2])/(1-a.mean))*sum.a
  
  piie.sp.ind <- data.y - psi.sp.ind
  
  piie.sp <- mean(data.y) - mean(psi.sp.ind)
  
  score.sp <- cbind( model.matrix.a*c(data.a - a.mean),
                     model.matrix.y*c(data.y - y.mean),
                     (piie.sp.ind - piie.est))
  
  deriv <- matrix(1,1,n)%*%score.sp
  
  return(deriv)
}



piie.nie.variance.function.cont <- function(cov.vals.all,exposure.data,theta.hat,beta.hat,alpha.hat,astar,interaction){
  
  n <- length(exposure.data)
  
  mean.cov.all <- apply(cov.vals.all,2,mean)
  mean.exposure <- mean(exposure.data)
  mean.cov.exposure <- apply(exposure.data*cov.vals.all,2,mean)
  
  
  if (interaction == 1){psi <- (theta.hat[1] + theta.hat[3]*beta.hat[1] + theta.hat[3]*beta.hat[2]*astar
                                + (theta.hat[2] + theta.hat[4]*beta.hat[1] + theta.hat[4]*beta.hat[2]*astar)*mean.exposure
                                + (theta.hat[3]*beta.hat[3:length(beta.hat)] + theta.hat[5:length(theta.hat)])%*%t(t(mean.cov.all))
                                + theta.hat[4]*beta.hat[3:length(beta.hat)]%*%t(t(mean.cov.exposure)))
  } else {psi <- (theta.hat[1] + theta.hat[3]*beta.hat[1] + theta.hat[3]*beta.hat[2]*astar
                  + theta.hat[2]*mean.exposure
                  + (theta.hat[3]*beta.hat[3:length(beta.hat)] + theta.hat[4:length(theta.hat)])%*%t(t(mean.cov.all)))}
  
  if(interaction == 1){est.piie <- theta.hat[3]*beta.hat[2]*mean.exposure + theta.hat[4]*beta.hat[2]*mean.exposure
  } else {est.piie <- theta.hat[3]*beta.hat[2]*mean.exposure   }
  
  if(interaction ==1){var.piie <- ( (var(exposure.data)/n)*(beta.hat[2]^2)*(theta.hat[3] + theta.hat[4])^2 
                                    + (mean.exposure^2)*((theta.hat[3] + theta.hat[4])^2)*vcov(fit.z)[2,2]
                                    + mean.exposure*beta.hat[2]*(mean.exposure*beta.hat[2]*vcov(fit.y)[3,3] + mean.exposure*beta.hat[2]*vcov(fit.y)[3,4])
                                    + mean.exposure*beta.hat[2]*(mean.exposure*beta.hat[2]*vcov(fit.y)[3,4] + mean.exposure*beta.hat[2]*vcov(fit.y)[4,4]) )
  } else { var.piie <- ( ((mean.exposure*theta.hat[3])^2)*vcov(fit.z)[2,2]
                         + ((mean.exposure*beta.hat[2])^2)*vcov(fit.y)[3,3]
                         + ((beta.hat[2]*theta.hat[3])^2)*var(exposure.data)/n )  }
  
  output <- cbind(psi,est.piie,var.piie)
  colnames(output) <- c("PSI","PIIE","Var PIIE")
  
  return(output)
  
}



